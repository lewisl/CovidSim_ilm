
###################################################
# data structures for vaccination
###################################################

@Base.kwdef mutable struct Vaccineparams
    name::Symbol 
    reqdshots::Int
    delay2ndshot::Union{Int, Nothing}   # days until 2nd shot (probability less important)
    halflife::Int  # days to 50% decline in effectiveness
    sendrisk::Vector{Float64}
    recvrisk::Vector{Float64}
    recvrisk_reduction::Dict{Symbol, Float64}
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        Vaccineparams(vd::Dict)=(
            Vaccineparams(
                name                     = Symbol(vd[:name]),
                reqdshots                = vd[:reqdshots],
                delay2ndshot             = vd[:delay2ndshot],
                halflife                 = vd[:halflife],
                sendrisk                 = vd[:sendrisk],
                recvrisk                 = vd[:recvrisk],
                recvrisk_reduction       = vd[:recvrisk_reduction], 
                )
        )


"""
    mutable struct Vaxinclude

Describes one vaccine included in a vaccination schedule. A dict holds instances of
this struct for each vaccine type included in a vaccination schedule.
"""
@Base.kwdef mutable struct Vaxinclude
    mix::Float64
    doses::Int
    pct2ndshot::Float64
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        function Vaxinclude(vi::Dict)
            mix = vi[:mix]
            doses = vi[:doses]
            pct2ndshot=vi[:pct2ndshot]

            Vaxinclude(mix=mix, doses=doses, pct2ndshot=pct2ndshot)
        end


@Base.kwdef mutable struct Vaxsched
    vaxesincluded::Dict{Symbol, Vaxinclude}
    dayrange::UnitRange{Int64}
    targetpct::Float64
    filterfunc::Function  # discretionary criteria for who is included in the vaccine schedule
    shotmode::Symbol   # one of :first, :second, :all, :booster.  best to use :all. others force only single shot
    pattern::Vector{Float64} # = [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]   
        # pattern defines 11 point on a piecewise linear "curve" of the pace of vaccination as the
        # percentage of a given population that receives vaccination-> must sum to 1.0
    pctfunc::Function

end
        # external method to accept inputs from a dict
        function Vaxsched(vs::Dict) 
            dayrange = vs[:dayrange][1]:vs[:dayrange][2]
            targetpct = vs[:targetpct]
            pattern = vs[:pattern]
            shotmode = Symbol(vs[:shotmode])

            vx =Vaxsched(
                vaxesincluded = Dict(k => Vaxinclude(v) for (k,v) in vs[:vaxesincluded]),
                dayrange =  dayrange,
                targetpct = targetpct,
                filterfunc = basevaxfilterfunc,   # excluded to enable later definition
                shotmode = shotmode,
                pattern = pattern,
                pctfunc = makevaxpctperdayfn(dayrange, targetpct, pattern; shotmode=shotmode) # returns a function that takes a day as input
                )   

            return vx
        end

        # internal method to verify mix, use defaults for filterfunc, pattern, shotmode and pctfunc
        function Vaxsched(;
            vaxesincluded::Dict{Symbol, Vaxinclude},
            dayrange, 
            targetpct, 
            filterfunc=basevaxfilterfunc,
            pattern=[0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01], 
            shotmode=:all,
            pctfunc=makevaxpctperdayfn(dayrange, targetpct, pattern; shotmode=:all)
            )

            vaxmix = [v.mix for v in values(vaxesincluded)]
            @assert sum(vaxmix) == 1.0 "Sum of mixes for vaccines does not equal 1: $vaxmix"

            new(vaxesincluded,
                dayrange,
                targetpct,
                filterfunc,
                shotmode,
                pattern,
                pctfunc
                )
        end



###########################################################
# setup code for vaccination
###########################################################


# TODO: a container for multiple concurrent vax schedules
# TODO: setup to build both vaccines and vax schedules


# this works
function build_vaxset(vaccinefilename; paramdir="../parameters")

    vaccinefiles = YAML.load_file(joinpath(paramdir, vaccinefilename); dicttype=Dict{Symbol,Any})

    vaxset = Dict()
    for vax in keys(vaccinefiles)
        vparamsdict = YAML.load_file(joinpath(paramdir, "vaccine_parameters", vaccinefiles[vax][:directory_name],
            vaccinefiles[vax][:infect_fname]), dicttype=Dict{Symbol, Any})
        vparams = Vaccineparams(vparamsdict)
        vtrans = setup_dt(joinpath("../parameters", "vaccine_parameters", vaccinefiles[vax][:directory_name], 
                    vaccinefiles[vax][:transition_fname]))
        vaxset[vax] = Dict(:params => vparams, :transition => vtrans)
    end
    return vaxset
end


function build_vaxschedset(schedfiles; paramdir="../parameters", scheddir="vaccine_schedule")

    schedpath = joinpath(paramdir, scheddir)
    vaxschedset = Dict{Symbol, Vaxsched}()

    for schedfile in schedfiles
        schedname = first(splitext(schedfile))
        vaxscheddict = YAML.load_file(joinpath(schedpath, schedfile), dicttype=Dict{Symbol, Any})
        vaxschedset[Symbol(schedname)] = Vaxsched(vaxscheddict)
    end

    return vaxschedset
end



function build_vaxschedset(; paramdir="../parameters", scheddir="vaccine_schedule")
    fnames = readdir(joinpath(paramdir, scheddir), join=true)
    fnames = filter(isfile, fnames)
    fnames = filter(fn -> (splitext(fn)[2] == ".yml"), fnames)
    fnames = map(fn->last(splitpath(fn)), fnames)  # get rid of the path--keep only filename

    build_vaxschedset(fnames, paramdir=paramdir, scheddir=scheddir)
end


"""
    makevaxfn(dayrange, pattern)

    returns: function pctperday(day)

Create a function that implements the vaccination schedule for a given vaccine type.
This function, when called with a simulation day, returns the percentage of the
target population to receive the vaccine on that day.
"""
function makevaxpctperdayfn(dayrange, targetpct, pattern; shotmode=:all)
    schedlength = length(dayrange)
    distribscale = pattern .* ((length(pattern)-1)/schedlength)
    interp = LinearInterpolation(0:length(distribscale)-1, distribscale)
    startday = dayrange.start; endday1 = dayrange.stop
    function vaxpctperday(day)
        p = (day - startday + 1) * round((length(interp)-1) / schedlength, digits=4)
        p = p < length(pattern) - 1 ? p : length(pattern) - 1 
        return interp(p) * targetpct
    end
end


function vaccinate!(locdat, vxschedset, contactable_idx, vaxset)
        
    today = day_ctr[:day]

    for vxsched in values(vxschedset)

        # lots of setup before we process the vaccine schedule vxsched
        vaxesincluded = collect(keys(vxsched.vaxesincluded))

        dayrange = vxsched.dayrange
        delay2ndshot = Dict(v=>vaxset[v][:params].delay2ndshot for v in vaxesincluded)
        dvals = filter(!isnothing, collect(values(delay2ndshot)))
        maxdelay = length(dvals) == 0 ? 0 : maximum(dvals)
        stop =  maxdelay == 0 ? dayrange.stop : dayrange.stop + maxdelay

        # shortcircuit the whole shebang for this schedule 
        if (today < dayrange.start) | (today > stop)
            continue
        end

        # schedule or vaccine parameters
        filterfunc = vxsched.filterfunc
        shotmode = vxsched.shotmode # values in :first, :second, :all, :booster
        pctfunc = vxsched.pctfunc
        vaxprops = vxsched.vaxesincluded
        reqdshots = Dict(v => vaxset[v][:params].reqdshots for v in vaxesincluded)
        pct2ndshot = Dict(k => v.pct2ndshot for (k,v) in vxsched.vaxesincluded)
        mix = [v.mix for v in values(vxsched.vaxesincluded)]
        doses = Dict(k => v.doses for (k,v) in vxsched.vaxesincluded)  # doses available
        
        # columns
        vaxstatuscol = locdat.vaxstatus
        vaxdaycol = locdat.vaxday
        vaxrcvdcol = locdat.vaxrcvd
        statuscol = locdat.status
        condcol = locdat.cond
        agegrpcol = locdat.agegrp

        # how many shots to give today?
        shotsremaining = floor(Int, vxsched.pctfunc(today) * length(locdat))   # pct times population


        # TODO we need a way in the vax cases to kill a schedule: reset dayrange to 1:1
        # TODO need to handle booster shots

        # loop across people who are not dead
        for p in contactable_idx  

            # are there still shots for this schedule?
            shotsremaining < 1 && break  

            # short circuit out of person loop if no more doses left of any vaccine 
            sum(values(doses)) <= 0 && break

            if vaxstatuscol[p] == :full      
                # DO NOTHING: this person gets no more shots
                
            elseif vaxstatuscol[p] == :none  # maybe give the first shot

                # which vaccine to give?
                vxnum = categorical_sim(mix)  # our first choice, if available
                
                vaxchoice = :none
                for vx in union([vaxesincluded[vxnum]], vaxesincluded) # puts vxnum first
                    if doses[vx] > 0  # if never > 0, exhaust the loop, and vax choice remains :none
                        vaxchoice = vx
                        break
                    end
                end
                if vaxchoice != :none   # we found doses to give
                    shotsremaining -= 1
                    doses[vaxchoice] -= 1
                    # update person's traits
                    vaxrcvdcol[p] = [vaxchoice]
                    vaxdaycol[p] = [today]    
                    if reqdshots[vaxchoice] > 1
                        vaxstatuscol[p] = :first
                    else
                        vaxstatuscol[p] = :full
                    end
                end

            elseif (vaxstatuscol[p] == :first) | (vaxstatuscol[p] == :multiple)
                vaxchoice = last(vaxrcvdcol[p]) # assume we don't mix vaccines for multiple shots

                if reqdshots[vaxchoice] > 1
                    if doses[vaxchoice] > 0
                        # is it time for the next shot?
                        prev_date = last(vaxdaycol[p])
                        if (today - prev_date) >= delay2ndshot[vaxchoice]
                            # will this person get the 2nd shot?
                            donext = Bool(binomial_one_sample(1,pct2ndshot[vaxchoice]))
                            if donext
                                shotsremaining -= 1
                                doses[vaxchoice] -= 1
                                # update person's traits
                                push!(vaxrcvdcol[p], vaxchoice)
                                push!(vaxdaycol[p], today)
                                if length(vaxrcvdcol[p])  >= reqdshots[vaxchoice]
                                    vaxstatuscol[p] = :full
                                else
                                    vaxstatuscol[p] = :multiple
                                end
                            end
                        end
                    end
                end
            else   
                # more than one shot received, but not :full--not sure how to use
                @assert false "vax status conditions failed: no condition satisfied"
            end

        end  # for p
    end  # for vxsched

end


function basevaxfilterfunc(p, agegrpcol)
    (agegrpcol[p] != age0_19)
end