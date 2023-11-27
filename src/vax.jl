######################################################################
#  vax.jl: data structures, setup, and implementation for vaccination
######################################################################


# TODO
#   add a way to filter who gets vaccinated in the vaxschedule (or a case or in the call to vaccinate?)


###################################################
# data structures for vaccination
###################################################

@Base.kwdef struct Vaccineparams  
    reqdshots::Int
    delay2ndshot::Union{Int, Nothing}   # days until 2nd shot (probability less important)
    delaybooster::Union{Int, Nothing}
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
                delaybooster             = vd[:delaybooster],
                halflife                 = vd[:halflife],
                effectiveness            = vd[:effectiveness], 
                full_effect_days         = vd[:full_effect_days],
                day1_effect              = vd[:day1_effect],
                infectfactor             = Dict(k => convert(Float64, v) for (k,v) in vd[:infectfactor])
                )
        )
#       TODO write a Base.show() for this struct
#=
# this is used to handle a call to `print`
Base.show(io::IO, x::MyString) = print(io, x.s)
# this is used to show values in the REPL and when using IJulia
Base.show(io::IO, m::MIME"text/plain", x::MyString) = print(io, x.s)
=#

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
    pctboost::Float64
    alternate::Array{String}
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        function Vaxinclude(vi::Dict)
            mix             = vi[:mix]
            starting_doses  = vi[:starting_doses]
            pct2ndshot      = vi[:pct2ndshot]
            pctboost        = vi[:pctboost]
            alternate       = vi[:alternate]

            Vaxinclude(mix=mix, starting_doses=starting_doses, pct2ndshot=pct2ndshot, pctboost=pctboost, alternate=alternate)
        end


@Base.kwdef struct Vaxsched  
    vaxesincluded::Dict{Symbol, Vaxinclude}
    dayrange::UnitRange{Int64}
    targetpct::Float64
    filtervec::Vector  # discretionary criteria for who is included in the vaccine schedule
    shotmode::Symbol   # one of :first, :second, :all, :booster.  best to use :all. others force only single shot
    pattern::Vector{Float64} # = [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]   
        # pattern defines 11 point on a piecewise linear "curve" of the pace of vaccination as the
        # percentage of a given population that receives vaccination-> must sum to 1.0
    spreadfunc::Function
end

        # external method to accept inputs from a dict
        function Vaxsched(vs_dict::Dict) 
            vaxesincluded = Dict(k => Vaxinclude(v) for (k, v) in vs_dict[:vaxesincluded])
            dayrange = vs_dict[:dayrange][1]:vs_dict[:dayrange][2]
            targetpct = vs_dict[:targetpct]
            filtervec = symbol2agegrp.(vs_dict[:filtervec])
            shotmode = Symbol(vs_dict[:shotmode])
            pattern = vs_dict[:pattern]
            
            vaxmix = [v.mix for v in values(vaxesincluded)]
            @assert sum(vaxmix) == 1.0 "Sum of mix for vaccines does not equal 1: $vaxmix"
            @assert 0.0 <= targetpct <= 1.0 "targetpct must be in [0.0, 1.0], got: $targetpct"

            Vaxsched(
                vaxesincluded = vaxesincluded,
                dayrange =  dayrange,
                targetpct = targetpct,
                filtervec = filtervec,   # excluded to enable later definition
                shotmode = shotmode,
                pattern = pattern,
                spreadfunc = genvaxspreadfunc(dayrange, targetpct, pattern; shotmode=shotmode) 
                )   
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
    vaxlist = collect(keys(vaccines))

    for vax in vaxlist
        vaxset[vax] = Vaccineparams(vaccines[vax])  
    end

    return vaxset, vaxlist
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
    
    vaxscheds = change_key_type(vaxscheds, f=Symbol)

    vaxschedset = Dict{Symbol, Vaxsched}()

    for sched in keys(vaxscheds)
        vaxschedset[Symbol(sched)] = Vaxsched(vaxscheds[sched])
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


###########################################################
# run vaccination
###########################################################


"""
Give people shots!

Determine the number of doses for each vaccine available to administer today.
Call doshots! to select recipients for each dose.
"""
@inline function vaccinate!(locdat, vaxschedset, vaxset, whichvaxscheds)
    
    vaxscheds = setvaxscheds(whichvaxscheds, vaxschedset)  # return vector of symbols or nothing

    isnothing(vaxscheds) && return nothing  # exit the function 
        
    today = DAY_CTR[:day]

    @inbounds @fastmath for schedname in Symbol.(vaxscheds)      
        vxsched = vaxschedset[schedname]   # vaxshedset is Dict{Symbol, Vaxsched} where Symbol is Symbol(schedname)

        # setup this schedule
        vaxprops = vxsched.vaxesincluded  # a dict of vaccine symbol to the struct Vaxinclude
                 # that contains mix, doses, starting_doses, pct2ndshot, pctboost, alternate, booster
        vaxesincluded = collect(keys(vaxprops))

        dayrange = vxsched.dayrange
        delay2ndshot = Dict(v=>vaxset[v].delay2ndshot for v in vaxesincluded) # per vax in the schedule
        delaybooster = Dict(v=>vaxset[v].delaybooster for v in vaxesincluded)
        maxdelay = mapreduce(v->vaxset[v].delay2ndshot, max, vaxesincluded)
        stop =  dayrange.stop + maxdelay

        # shortcircuit the whole shebang for this schedule 
        if (today < dayrange.start) | (today > stop)
            continue
        end

        # schedule parameters
        filtervec    = vxsched.filtervec  # contains allowed agegrps
        shotmode      = vxsched.shotmode      # values in :first, :second, :all, :booster   TODO we are not using this yet
        spreadfunc    = vxsched.spreadfunc
        pct2ndshot    = Dict(k => v.pct2ndshot for (k,v) in vaxprops)   # per vax
        pctboost      = Dict(k => v.pctboost for (k,v) in vaxprops)   # per vax
        mix           = [v.mix for v in values(vaxprops)]

        # vax parameters
        reqdshots = Dict(v => vaxset[v].reqdshots for v in vaxesincluded)
        starting_doses = Dict(vi => vaxprops[vi].starting_doses for vi in vaxesincluded)  
        doses_today = Dict(vi => floor(Int, spreadfunc(today) * starting_doses[vi]) for vi in vaxesincluded)   # pct times accessible population

        # who is eligible to receive a vaccine?
            c_status = locdat.status
            c_recovday = locdat.recovday
            c_agegrp = locdat.agegrp

            filt_unexp = [x === :unexposed for x in c_status]  # comprehension fastest compared to map, loop, or . syntax
            # must use short-circuit && for and === :recovered->guarantees that element of c_recovday is not empty
            filt_recov = [(c_status[i] === :recovered) && (last(c_recovday[i]) < (today - 14)) for i in eachindex(c_status)]
            filt_agegrp = [in(x, filtervec) for x in c_agegrp]
            
            vaxable_idx = findall(filt_agegrp .&& (filt_unexp .| filt_recov))  # must be in agegrps and either unexposed or recovered
        
        doshots!(locdat,                        
                  vaxprops, vaxesincluded, reqdshots, pct2ndshot, pctboost, mix, delay2ndshot, delaybooster,   
                  vaxable_idx, doses_today, filtervec, today) 

    end  
end


"""
For input of :all, :none or a single symbol return a vector of vaxscheds as symbols.
"""
@inline function setvaxscheds(whichvaxsched::Symbol, vaxschedset)
    if whichvaxsched === :all
        return collect(keys(vaxschedset))
    elseif vaxscheds === :none
        return nothing
    else  # a single symbol turned into an array
        return [whichvaxsched]
    end
end

@inline function setvaxscheds(whichvaxsched::Vector{Symbol}, vaxschedset)
    return whichvaxsched
end


"""
Determine who gets a shot today and administer it; update population data.
"""
@inline function doshots!(locdat,                           # arrays to update
                  vaxprops, vaxesincluded, reqdshots, pct2ndshot, pctboost, mix, delay2ndshot, delaybooster,  # vaccine characteristics
                  vaxable_idx, doses_today, filtervec, today)      # agegrp                         # people and simulation today

    avail_doses = mapreduce(vi->doses_today[vi], +, keys(vaxprops))
    avail_doses <= 0  && return

    @inbounds @fastmath for p in shuffle!(vaxable_idx)

        # break out if no more doses left of any vaccine 
        # avail_doses = mapreduce(vi->doses_today[vi], +, keys(vaxprops))
        avail_doses <= 0 && break  

        person = locdat[p]
    
        if person.vaxstatus === :none  # maybe give the first shot

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
                avail_doses -= 1

                # update person's traits
                person.vaxrcvd = [vaxchoice]
                person.vaxday = [today]       # because this is first shot
                if reqdshots[vaxchoice] > 1
                    person.vaxstatus = :first
                else
                    person.vaxstatus = :full
                end
            end
                            # TODO add separate parameter for delaybooster
        elseif person.vaxstatus === :first
            vaxchoice = last(person.vaxrcvd) # assume we try not to mix vaccines for multiple shots

            if vaxprops[vaxchoice].doses > 0  # we have this vaccine in our remaining daily allotment
                # is it time for the next shot?
                prev_date = last(person.vaxday)
                if (today - prev_date) >= delay2ndshot[vaxchoice]
                    # will this person get another shot?  (based on pct2ndshot parameter input)
                    dotwo = Bool(binomial_one_sample(1, pct2ndshot[vaxchoice]))
                    if dotwo
                        vaxprops[vaxchoice].doses -= 1    # reduce the supply
                        doses_today[vaxchoice] -= 1       # reduce today's allotment
                        avail_doses -= 1

                        # update person's traits
                        push!(person.vaxrcvd, vaxchoice)
                        push!(person.vaxday, today)
                        person.vaxstatus = :full
                    end  # if dotwo
                end  # time for next shot
            end  # doses > 0

        elseif (person.vaxstatus === :booster) | (person.vaxstatus === :full)
            vaxchoice = last(person.vaxrcvd) # assume we try not to mix vaccines for multiple shots

            if vaxprops[vaxchoice].doses > 0  # we have this vaccine in our remaining daily allotment
                # is it time for the next shot?
                prev_date = last(person.vaxday)
                if (today - prev_date) >= delaybooster[vaxchoice]
                    # will this person get another shot?  (based on pct2ndshot parameter input)
                    domore = Bool(binomial_one_sample(1, pctboost[vaxchoice]))
                    if domore
                        vaxprops[vaxchoice].doses -= 1    # reduce the supply
                        doses_today[vaxchoice] -= 1       # reduce today's allotment
                        avail_doses -= 1

                        # update person's traits
                        push!(person.vaxrcvd, vaxchoice)
                        push!(person.vaxday, today)
                        person.vaxstatus = :booster     # assume everything over full is a booster...
                    end  # if domore
                end  # time for next shot
            end  # doses > 0

        end  # if c_vaxstatus -> time for another shot after the 1st shot

    end  # for p
end

