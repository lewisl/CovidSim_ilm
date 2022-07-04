# TODO
#   add a way to filter who gets vaccinated in the vaxschedule (or a case or in the call to vaccinate?)


###################################################
# data structures for vaccination
###################################################

const vaxlist = Symbol[]  # filled as vaxset is built

@Base.kwdef struct Vaccineparams  # mutable  ??
    reqdshots::Int
    delay2ndshot::Union{Int, Nothing}   # days until 2nd shot (probability less important)
    halflife::Int  # days to 50% decline in effect
    effectiveness::Dict{Symbol, Dict{Symbol, Float64}}
    full_effect_days::Int
    day1_effect::Float64
    infectfactor::Dict{Symbol, Float64}
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        Vaccineparams(vd::Dict)=(
            Vaccineparams(
                reqdshots                = vd[:reqdshots],
                delay2ndshot             = vd[:delay2ndshot],
                halflife                 = vd[:halflife],
                effectiveness            = vd[:effectiveness], 
                full_effect_days         = vd[:full_effect_days],
                day1_effect              = vd[:day1_effect],
                infectfactor             = Dict(k => convert(Float64, v) for (k,v) in vd[:infectfactor])
                )
        )


"""
    struct Vaxinclude

Describes one vaccine included in a vaccination schedule. A dict holds instances of
this struct for each vaccine type included in a vaccination schedule.
"""
@Base.kwdef mutable struct Vaxinclude   
    mix::Float64
    starting_doses::Int
    doses::Int=0   # set equal to starting doses when initializing beginning of simulation run
    pct2ndshot::Float64
    alternate::Array{String}
    booster::Bool
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        function Vaxinclude(vi::Dict)
            mix             = vi[:mix]
            starting_doses  = vi[:starting_doses]
            pct2ndshot      = vi[:pct2ndshot]
            alternate       = vi[:alternate]
            booster         = vi[:booster]

            Vaxinclude(mix=mix, starting_doses=starting_doses, pct2ndshot=pct2ndshot, alternate=alternate, booster=booster)
        end


@Base.kwdef struct Vaxsched  
    vaxesincluded::Dict{Symbol, Vaxinclude}
    dayrange::UnitRange{Int64}
    targetpct::Float64
    filterfunc::Vector  # discretionary criteria for who is included in the vaccine schedule
    shotmode::Symbol   # one of :first, :second, :all, :booster.  best to use :all. others force only single shot
    pattern::Vector{Float64} # = [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]   
        # pattern defines 11 point on a piecewise linear "curve" of the pace of vaccination as the
        # percentage of a given population that receives vaccination-> must sum to 1.0
    spreadfunc::Function
end

        # external method to accept inputs from a dict
        function Vaxsched(vs::Dict) 
            dayrange = vs[:dayrange][1]:vs[:dayrange][2]
            targetpct = vs[:targetpct]
            pattern = vs[:pattern]
            filterfunc = symbol2agegrp.(vs[:filterfunc])
            shotmode = Symbol(vs[:shotmode])

            vaxesincluded = Dict(k => Vaxinclude(v) for (k,v) in vs[:vaxesincluded])
            vaxmix = [v.mix for v in values(vaxesincluded)]

            @assert sum(vaxmix) == 1.0 "Sum of mix for vaccines does not equal 1: $vaxmix"
            @assert 0.0 <= targetpct <= 1.0 "targetpct must be in [0.0, 1.0], got: $targetpct"

            vx =Vaxsched(
                vaxesincluded = vaxesincluded,
                dayrange =  dayrange,
                targetpct = targetpct,
                filterfunc = filterfunc,   # excluded to enable later definition
                shotmode = shotmode,
                pattern = pattern,
                spreadfunc = genvaxspreadfunc(dayrange, targetpct, pattern; shotmode=shotmode) 
                )   

            return vx
        end


###########################################################
# setup code for vaccination
###########################################################


function build_vaxset(vaccinefilename, paramdir)
    vaccines = YAML.load_file(joinpath(paramdir, vaccinefilename); dicttype=Dict{Symbol,Any})
    build_vaxset(vaccines)
end


function build_vaxset(vaccines)

    vaxset = Dict{Symbol, Vaccineparams}()
    for vax in keys(vaccines)
        vaxset[vax] = Vaccineparams(vaccines[vax])  
    end

    # clean vaxlist
    l = length(vaxlist)
    if l > 0
        deleteat!(vaxlist, collect(1:l))
    end

    if isempty(vaxlist)
        for k in keys(vaxset)
            push!(vaxlist, k)   # this is a module global variable. Forgive me for I have sinned--except it makes sense...
        end
    end

    return vaxset
end


function build_vaxschedset(scheddir, paramdir)
    fnames = readdir(joinpath(paramdir, scheddir), join=true)
    fnames = filter(isfile, fnames)
    fnames = filter(fn -> (splitext(fn)[2] == ".yml"), fnames)
    fnames = map(fn->last(splitpath(fn)), fnames)  # get rid of the path--keep only filename

    # build_vaxschedset(fnames, paramdir=paramdir, scheddir=scheddir)
    schedpath = joinpath(paramdir, scheddir)
    vaxschedset = Dict{Symbol, Vaxsched}()

    for schedfile in fnames
        schedname = first(splitext(schedfile))
        vaxscheddict = YAML.load_file(joinpath(schedpath, schedfile), dicttype=Dict{Symbol, Any})
        vaxschedset[Symbol(schedname)] = Vaxsched(vaxscheddict)
    end

    return vaxschedset
end


function build_vaxschedset(vaxscheds)  # input is a Dict{Any, Any}
    
    vaxschedset = Dict{Symbol, Vaxsched}()

    for sched in keys(vaxscheds)
        vaxscheddict = YAML.load(vaxscheds[sched], dicttype=Dict{Symbol, Any})
        vaxschedset[Symbol(sched)] = Vaxsched(vaxscheddict)
    end

    return vaxschedset
end

"""
    genvaxspreadfunc(dayrange, targetpct, pattern; shotmode=:all)

Create a function that implements the vaccination schedule for a given vaccine type.
This function, when called with a simulation day, returns the percentage of the
target population to receive the vaccine on that day.
"""
function genvaxspreadfunc(dayrange, targetpct, pattern; shotmode=:all)
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


"""
Give people shots!
"""
@inline function vaccinate!(locdat, vaxschedset, vaxset, vaxscheds::Vector{Symbol})
        
    today = day_ctr[:day]

    for schedname in Symbol.(vaxscheds)      
        vxsched = vaxschedset[schedname]   # vaxshedset is Dict{Symbol, Vaxsched} where Symbol is Symbol(schedname)

        # setup this schedule
        vaxprops = vxsched.vaxesincluded  # a dict of vaccine symbol to the struct Vaxinclude
                 # that contains mix, doses, starting_doses, pct2ndshot, alternate, booster
        vaxesincluded = collect(keys(vaxprops))

        dayrange = vxsched.dayrange
        delay2ndshot = Dict(v=>vaxset[v].delay2ndshot for v in vaxesincluded) # per vax in the schedule
        maxdelay = mapreduce(v->vaxset[v].delay2ndshot, max, vaxesincluded)
        stop =  dayrange.stop + maxdelay

        # shortcircuit the whole shebang for this schedule 
        if (today < dayrange.start) | (today > stop)
            continue
        end

        # schedule parameters
        filterfunc    = vxsched.filterfunc  # contains allowed agegrps
        shotmode      = vxsched.shotmode      # values in :first, :second, :all, :booster   TODO we are not using this yet
        spreadfunc    = vxsched.spreadfunc
        pct2ndshot    = Dict(k => v.pct2ndshot for (k,v) in vaxprops)   # per vax
        mix           = [v.mix for v in values(vaxprops)]

        # vax parameters
        reqdshots = Dict(v => vaxset[v].reqdshots for v in vaxesincluded)
        
        # people columns
        vaxstatuscol    = locdat.vaxstatus
        vaxdaycol       = locdat.vaxday
        vaxrcvdcol      = locdat.vaxrcvd
        agegrpcol       = locdat.agegrp

        
        starting_doses = Dict(vi => vaxprops[vi].starting_doses for vi in vaxesincluded)  

        doses_today = Dict(vi => floor(Int, spreadfunc(today) * starting_doses[vi]) for vi in vaxesincluded)   # pct times accessible population

        vaxable_idx = findall((locdat.status .== unexposed) 
                                .| ((locdat.status .== recovered) .& (last.(locdat.recovday) .< today - 14))
                                .& (in.(locdat.agegrp, [filterfunc]))   # horrible syntax! (for included age groups)
                                ) 
        
        doshots!(vaxrcvdcol, vaxdaycol, vaxstatuscol,  
                  vaxprops, vaxesincluded, reqdshots, pct2ndshot, mix, delay2ndshot,    
                  vaxable_idx, doses_today, agegrpcol, filterfunc, today)

    end  # for schedname
end


"""
Front-end method for vaccinate!Allows vaxscheds argument to be :all, :none or a single symbol for a specific vaxsched.
"""
@inline function vaccinate!(locdat, vaxschedset, vaxset, vaxscheds::Symbol)
    if vaxscheds === :all
        vaxscheds = collect(keys(vaxschedset))
        vaccinate!(locdat, vaxschedset, vaxset, vaxscheds)

    elseif vaxscheds === :none
        # don't do anything
    else  # a single symbol turned into an array
        vaccinate!(locdat, vaxschedset, vaxset, [vaxscheds])  
    end
end



@inline function doshots!(vaxrcvdcol, vaxdaycol, vaxstatuscol,         # arrays to update
                  vaxprops, vaxesincluded, reqdshots, pct2ndshot, mix, delay2ndshot,  # vaccine characteristics
                  vaxable_idx, doses_today, agegrpcol, filterfunc, today)                               # people and simulation today

    # break out if all the people in this schedule today have been fully vaccinated
    # lreturn  nothing      # sum(values(doses_today)) < 1 && break  

    for p in shuffle(vaxable_idx)

        # break out if no more doses left of any vaccine 
        avail_doses = mapreduce(vi->doses_today[vi], +, keys(vaxprops))
        avail_doses <= 0 && break  
    
        if vaxstatuscol[p] == :none  # maybe give the first shot

            # which vaccine to give?
            vxnum = categorical_sim(mix)  # our first choice, if available
            vaxchoice = vaxesincluded[vxnum]

            if vaxprops[vaxchoice].doses < 1  # out of first choice!
                for alt in vaxprops[vaxchoice].alternate
                    alt = Symbol(alt)
                    if vaxprops[alt].doses > 0
                        vaxchoice = alt
                        break
                    else
                        vaxchoice=:none
                    end
                end
            end

            if vaxchoice != :none   # we found doses to give
                vaxprops[vaxchoice].doses -= 1   # reduce the supply
                doses_today[vaxchoice] -= 1      # reduce today's allotment

                # update person's traits
                vaxrcvdcol[p] = [vaxchoice]
                vaxdaycol[p] = [today]       # because this is first shot
                if reqdshots[vaxchoice] > 1
                    vaxstatuscol[p] = :first
                else
                    vaxstatuscol[p] = :full
                end
            end
                            # TODO separate logic branch for boosters using NEW booster pct input
        elseif (vaxstatuscol[p] === :first) | (vaxstatuscol[p] === :booster) | (vaxstatuscol[p] === :full)
            vaxchoice = last(vaxrcvdcol[p]) # assume we try not to mix vaccines for multiple shots

            if vaxprops[vaxchoice].doses > 0  # we have this vaccine in our remaining daily allotment
                # is it time for the next shot?
                prev_date = last(vaxdaycol[p])
                if (today - prev_date) >= delay2ndshot[vaxchoice]
                    # will this person get another shot?  (based on pct2ndshot parameter input)
                    dotwo = Bool(binomial_one_sample(1, pct2ndshot[vaxchoice]))
                    if dotwo
                        vaxprops[vaxchoice].doses -= 1    # reduce the supply
                        doses_today[vaxchoice] -= 1       # reduce today's allotment

                        # update person's traits
                        push!(vaxrcvdcol[p], vaxchoice)
                        push!(vaxdaycol[p], today)
                        if length(vaxrcvdcol[p])  == reqdshots[vaxchoice]
                            vaxstatuscol[p] = :full
                        elseif length(vaxrcvdcol[p])  >= reqdshots[vaxchoice]
                            vaxstatuscol[p] = :booster     # assume everything over full is a booster...
                        end
                    end  # if dotwo
                end  # time for next shot
            end  # doses > 0
        end  # if vaxstatuscol -> time for another shot after the 1st shot

    end  # for p
end


function basevaxfilterfunc(p, agegrpcol)
    (agegrpcol[p] != age0_19)
end