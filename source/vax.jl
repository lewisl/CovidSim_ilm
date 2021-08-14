
###################################################
# data structures for vaccination
###################################################

@Base.kwdef mutable struct Vaccineparams
    name::Symbol 
    shots::Int
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



@Base.kwdef mutable struct Vaxsched
    vaxinclude::Dict{Symbol, NamedTuple{(:mix, :doses, :pct2ndshot), Tuple{Float64, Int64, Float64}}}
    dayrange::UnitRange{Int64}
    targetpct::Float64
    filterfunc::Function  # discretionary criteria for who is included in the vaccine schedule
    shotmode::Symbol   # one of :first, :second, :all.  best to use :all. others force only single shot
    # pattern defines 11 point on a piecewise linear "curve" of the pace of vaccination as the
    # percentage of a given population that receives vaccination-> must sum to 1.0
    # [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]
    pattern::Vector{Float64} # = [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]         
    pctfunc::Function
end


        # this works
        """
        Method to build schedule of people to be vaccinated by percent by day.
        This will typically be done as part of a runcase rather than
        as part of initial setup. ???
        """
        function create_vaxsched(;vaxinclude::Dict{Symbol, NamedTuple{(:mix, :doses, :pct2ndshot), Tuple{Float64, Int64, Float64}}},
                                dayrange, targetpct, 
                                filterfunc=basevaxfilterfunc,
                                pattern=[0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01], 
                                shotmode=:all)

            vaxmix = [v.mix for v in values(vax)]
            @assert sum(vaxmix) == 1.0 "Sum of mixes for vaccines does not equal 1: $vaxmix"

            Vaxsched(
                vaxinclude=vaxinclude,
                dayrange=dayrange,
                targetpct=targetpct,
                filterfunc=filterfunc,
                shotmode=shotmode,
                pattern=pattern,
                pctfunc=makevaxpctperdayfn(dayrange, targetpct, pattern; shotmode=:all)
                )
        end


###########################################################
# setup code for vaccination
###########################################################


# TODO: a container for multiple concurrent vax schedules
# TODO: setup to build both vaccines and vax schedules
# TODO: a YAML file format for vax schedules


# this works
function setupvax(paramdir="../parameters")
    vaccinefiles = YAML.load_file(joinpath(paramdir, "vaccines.yml"); dicttype=Dict{Symbol,Any})

    pprint(vaccinefiles)

    vaxset = build_vaxset(vaccinefiles, paramdir)

    return vaxset
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

# this works
function build_vaxset(vaccinefiles, paramdir="../parameters")
    vaxset = Dict{Symbol, Dict{Symbol, Union{Vaccineparams, transitionT}}}()
    for vax in keys(vaccinefiles)
        vparams = YAML.load_file(joinpath(paramdir,vaccinefiles[vax][:directory_name],
            vaccinefiles[vax][:infect_fname]), dicttype=Dict{Symbol, Any})
        vparams = Vaccineparams(vparams)
        vtrans = setup_dt(joinpath("../parameters", vaccinefiles[vax][:directory_name], 
                    vaccinefiles[vax][:transition_fname]))
        vaxset[vax] = Dict()
        vaxset[vax][:params] = vparams
        vaxset[vax][:transition] = vtrans
    end
    return vaxset
end


function vaccinate!(locdat, vxschedset, contactable_idx, spreadset, dovax, dovariant)
        
    today = day_ctr[:day]

    for vxsched in vxschedset

        # shortcircuit the whole shebang for this vaccine in the schedule
        dayrange = vxsched.dayrange
        delay2ndshot = Dict(v=>spreadset[v].delay2ndshot for v in vaccines)
        maxstop2ndshot = dayrange.stop + maximum(values(delay2ndshot))

        if (today < dayrange.start) | (today > maxstop2ndshot)
            continue
        end

        # schedule parameters
        filterfunc = vxsched.filterfunc
        shotmode = vxsched.shotmode # values in :first, :second, :all
        pctfunc = vxsched.pctfunc
        vaxinclude = collect(keys(vxsched.vaxinclude))
        vaxprops = vxsched.vax
        reqdshots = Dict(v=>spreadset[v].reqdshots for v in vaccinesinsched)
        mix = [v.mix for v in values(vxsched.vax)]
        doses = Dict(k => v.doses for (k,v) in values.vax)  # doses available
        
        # columns
        vaxstatuscol = locdat.vaxstatus
        vaxdaycol = locdat.vaxday
        vaxrecdcol = locdat.vaxrecd
        statuscol = locdat.status
        condcol = locdat.cond
        agegrpcol = locdat.agegrp

        # how many shots to give today?
        shotsremaining = vxsched.pctfunc(today)


        # TODO we should stop if no more doses in the schedule
        # TODO we need a way in the vax cases to kill a schedule

        while shotsremaining > 0
            for p in contactable_idx  # people who are not dead

                if vaxstatuscol[p] == :full      
                    # DO NOTHING: this person gets no more shots
                elseif vaxstatuscol[p] == :none  # the first shot

                    # which vaccine to give?
                    vxnum = categorical_sim(mix)  # our first choice, if available
                    
                    vaxchoice = :none
                    for vx in union(vaxinclude[vxnum], vaxinclude) # put vxnum first
                        if doses[vx] > 0
                            vaxchoice = vx
                            break
                        end
                    end
                    if vaxchoice != :none   # we found doses to give
                        shotsremaining -= 1
                        doses[vaxchoice] -= 1
                        # update person's traits
                        vacrcvdcol[p] = [vaxchoice]
                        vaxdaycol[p] = [day]    
                        vaxstatuscol[p] = :first
                    end

                elseif vaxstatuscol[p] == :first  
                    vaxchoice = last(vaxrcvdcol[p]) # assume we don't mix vaccines for multiple shots
                    if doses[vaxchoice] > 0
                        shotsremaining -= 1
                        doses[vaxchoice] -= 1
                        # update person's traits
                        push!(vaxrcvdcol[p], thisvax)
                        push!(vaxdaycol[p], day)
                        if reqdshots[vaxchoice] == length(vaxrcvdcol[p])
                            vaxstatuscol[p] = :full
                        end
                    end
                else   
                    # more than one shot received, but not :full--not sure how to use
                end

            end  # for p
        end  # while shotsremaining
    end  # for vxsched

end


function basevaxfilterfunc(p, agegrpcol)
    (agegrpcol[p] != age0_19)
end