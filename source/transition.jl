####################################################
# transition.jl for ilm model
#     change status of folks in simulation:
#           transition
#           travel
####################################################

    
"""
    transition!(locdat, infect_idx, dectree)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. Finally,  
they move to recovered or dead.

locdat must be a population table for a single locale.
"""
@inline function transition!(locdat, infect_idx, infectset, transitionset, vaxset, dovax, dovariant; vaxfn! = noop, trvec = zeros(6))
        
    # aliases for person attribute columns--deref the named tuple once
    c_sickday = locdat.sickday
    c_cond = locdat.cond
    c_agegrp = locdat.agegrp
    c_vaxstatus = locdat.vaxstatus
    c_status = locdat.status
    c_variant = locdat.variant 
    c_recovday = locdat.recovday 
    c_vaxrcvd = locdat.vaxrcvd
    c_vaxday = locdat.vaxday

    if dovax
        vaxfn! = vaxfn! == noop ? vaxtransitioneffect! : vaxfn!  # last branch new vaxfn! was passed in
    end

    # TODO test variant of each person
    # TODO based on variant use adjustment of :base or :base
    transarray = transitionset[:base].tree   

    for p in infect_idx  # p for infected person    
        p_sickday = c_sickday[p]
        p_cond = c_cond[p]
        p_agegrp = c_agegrp[p]  # agegroup of person p = agegrp column of locale data, row p 

        # if person's agegrp, sickday, and condition match a transition stage
        transvec = has(getfield(transarray, Symbol(p_agegrp)), p_sickday, p_cond, trvec) 

        riskadjustments = riskadjust(infectset, vaxset, p;
                        c_vaxstatus=c_vaxstatus, c_status=c_status, c_variant=c_variant, 
                        c_recovday=c_recovday, c_vaxrcvd=c_vaxrcvd, c_vaxday=c_vaxday)   # TODO make a function specific to transition

        vaxfn!(transvec, riskadjustments) #vaxfn! will be function noop or function vaxtransitioneffect

        dotransition!(locdat, p, p_cond, transvec) # perform transition logic and update population table

    end  
end


@inline function has(agetr, sickday::Int, p_cond::condition, trvec)::Union{Vector{Float64}, Nothing}
    for trdef in agetr
        if trdef.sickday == sickday
            trvec[:] = trdef.transition[map2cond(p_cond), :]
            if sum(trvec) > 0.0
                return trvec
            end
        end
    end
    return nothing
end


# TODO: need to do effect of immunity, vax, variant
# this will be set to the variable vaxfn!
@inline function vaxtransitioneffect!(transvec, adjustments)
    if transvec === nothing
        return
    end

    combinedfactor = adjustments.variant * min(adjustments.recov, adjustments.vax)
    riskfactor = squashfunc(combinedfactor)  

    for c in (sick, severe, dead)
        transvec[map2transition(c)] *= riskfactor
    end

    correction = 1.0 / sum(transvec)
    transvec[:] .*= correction

end

"""
    dotransition!(locdat, p, p_cond, trvec::Union{Vector{Float64}, Nothing})

Transition an infected person to a new condition or status if called
with a transition vector (trvec) or increment
the number of days the person has been sick.
"""
@inline function dotransition!(locdat, p, p_cond, trvec::Union{Vector{Float64}, Nothing})
   
    if isnothing(trvec)

        locdat.sickday[p] += 1  

    else
        choice = categorical_sim(trvec) # which outcome based on probability...?

        # debugging
        if choice == 0
            println("debugging function dotransition!")
            println(trvec)
            @assert false
        end

        tocond = map2transition(choice)  # next condition or status

        # if locdat.sickday[p] >= 25
        #     println("$(day_ctr[:day]): agegrp: $(locdat.agegrp[p]) sickday: $(locdat.sickday[p]) from cond: $p_cond to cond: $tocond")
        # end

        if tocond == dead  
            locdat.deadday[p] = day_ctr[:day]
            locdat.status[p] = dead  # change the status
            # locdat.cond[p] = uninfected # change the condition--> kept to know what cause of death was
        elseif tocond == recovered
            push!(locdat.recovday[p], day_ctr[:day])
            locdat.status[p] = recovered
            # locdat.cond[p] = uninfected   # TODO decide if this makes sense--using this to maintain a history of past infection
        else   
            locdat.cond[p] = tocond   # change the condition = degree of sickness
            locdat.sickday[p] += 1    # advance number of days person has been sick
        end    
    end
end

"""
    riskadjust(infectset, vaxset, locdat, target) 

Adjust the risk of getting infected and the effect on transitioning through stages of the infection based
on vaccination, recovery from previous infection and the variant of the infection contracted by an individual.

Returns (variant=variantfactor, recov=recovfactor, vax=vaxfactor)
"""
@inline function tr_riskadjust(infectset, vaxset, locdat, target; spr=0)
    # calculate based on recovery date and variant half-life of partial immunity and immunity strength
    # initially using exponential decay   TODO add parameter for exponential or sigmoid decay

    # @bp

    today = day_ctr[:day]
    oneshotfactor = 0.85   # TODO yet another parameter to put somewhere...!

    # target person characteristics
    vaxstatus = locdat.vaxstatus[target]
    status = locdat.status[target]
    
    # recovery effect, based on variant of person's infection
    if status == recovered

        variant = locdat.variant[target][end]
        days_post_recov = today - locdat.recovday[target][end] #recovday is a vector of days--get the last one

        if days_post_recov > 0   # TODO should be an assert: does this run day of or day after recovery?
            # get the max immunity
            immstrength = infectset[variant].immunestrength

            # get the declined value
            immhalflife = infectset[variant].immunehalflife
            immdecline = lindecay(days_post_recov, immhalflife, 0.05)
            recovfactor = 1.0 - (immdecline * immstrength)
        end
    else
        recovfactor = 1.0  # no immunity effect from recovery
    end

    # vaccine effect (rise time and decay)
    if vaxstatus != :none
        # person's vaccine conditions
        vaxrcvd = locdat.vaxrcvd[target]
        vaxday = locdat.vaxday[target]
        @assert size(vaxrcvd, 1) == size(vaxday, 1) "Oh, no: lengths of vaxrcvd and vaxday not equal"
        days_post_vax = today - last(vaxday)

        # vaccine characteristics
        lastvax = last(vaxrcvd)
        halflife = vaxset[lastvax].halflife
        full_effect_days = vaxset[lastvax].full_effect_days
        infectfactor = vaxset[lastvax].infectfactor
     
        # based on variant: for transition--own variant; for spreading: spreader's variant
        variant = spr == 0 ? locdat.variant[target][end] : locdat.variant[spr][end]
        vaxeffect = vaxset[lastvax].effectiveness[vaxstatus][variant]

        # rise & decay
        vaxmod = vaxmodifier(full_effect_days, today, days_post_vax, halflife; rise_lower=0.5, decay_lower=0.05)
        vaxfactor = 1.0 - (vaxmod * vaxeffect * infectfactor)
    else
        vaxfactor = 1.0   # no immunity effect from vaccination
    end

    # variant effect (sender's variant affects infectiousness up or down)
    variantfactor = 1.0   # replace with actual calculation...


    return (variant=variantfactor, recov=recovfactor, vax=vaxfactor)

end


function map2transition(choice::Int) # faster than using a Dict, array, or tuple because few items

    if choice == 1  # most common
        recovered
    elseif choice == 2
        nil
    elseif choice == 3
        mild
    elseif choice == 4
        sick
    elseif choice == 5
        severe
    elseif choice == 6   # least common
        dead
    else
        @assert false "invalid condition or status integer value $choice"
    end
        
end

function map2transition(choice::Enum) # faster than using a Dict, array, or tuple because few items

    if choice == recovered  # most common
        1
    elseif choice == nil
        2
    elseif choice == mild
        3
    elseif choice == sick
        4
    elseif choice == severe
        5
    elseif choice == dead   # least common
        6
    else
        @assert false "invalid condition or status enum value $choice"
    end
        
end


function map2cond(x::Enum)
    if x==nil
        1
    elseif x==mild
        2
    elseif x==sick
        3
    elseif x==severe
        4
    else
        @assert false "invalid condition Enum $x"
    end
end

function map2cond(x::Int)
    if x==1
        nil
    elseif x==2
        mild
    elseif x==3
        sick
    elseif x==4
        severe
    else
        @assert false "invalid condition Int $x"
    end
end


"""
For a locale, randomly choose the number of people from each agegroup with
condition of {unexposed, infectious, recovered} who travel to each
other locale. Add to the travelq.
"""
function travelout!(fromloc, locales, rules=[])    # TODO THIS WON'T WORK ANY MORE!
    # 10.5 microseconds for 5 locales
    # choose distribution of people traveling by age and condition:
        # unexposed, infectious, recovered -> ignore sickday for now
    # TODO: more frequent travel to and from Major and Large cities
    # TODO: should the caller do the loop across locales?   YES
    travdests = collect(locales)
    deleteat!(travdests,findfirst(isequal(fromloc), travdests))
    bins = lim = length(travdests) + 1
    for agegrp in agegrps
        for cond in [unexposed, infectious, recovered]
            name = condnames[cond]
            for sickday in sickdays
                numfolks = sum(grab(cond, agegrp, sickday, fromloc)) # the from locale, all sickdays
                travcnt = floor(Int, gamma_prob(travprobs[agegrp]) * numfolks)  # interpret as fraction of people who will travel
                x = rand(travdests, travcnt)  # randomize across destinations
                bydest = bucket(x, vals=1:length(travdests))
                for dest in 1:length(bydest)
                    isempty(bydest) && continue
                    cnt = bydest[dest]
                    iszero(cnt) && continue
                    enqueue!(travelq, travitem(cnt, fromloc, dest, agegrp, sickday, name))
                end
            end
        end
    end
end


"""
Assuming a daily cycle, at the beginning of the day
process the queue of travelers from the end of the previous day.
Remove groups of travelers by agegrp, sickday, and condition
from where they departed.  Add them to their destination.
"""
function travelin!(dat=popdat)   # TODO THIS DOESN'T WORK ANYMORE
    while !isempty(travelq)
        g = dequeue!(travelq)
        cond = eval(Symbol(g.cond))
        minus!(g.cnt, cond, g.agegrp, g.sickday, g.from, dat=dat)
        plus!(g.cnt, cond, g.agegrp, g.sickday, g.to, dat=dat)
    end
end


"""
Return boolean filter of people currently in quarantine less
daily leakage, if applicable.
"""
function current_quar(locdat, leakage = .05)
    @assert 0.0 <= leakage <= 1.0 "leakage value must be between 0.0 and 1.0 inclusive"

    iq_filt = copy(locdat.quar)
    cnt = round(Int, sum(locdat.quar) * leakage)
    
    iq_filt[sample(findall(locdat.quar), cnt, replace=false)] .= false
    
    return iq_filt
end


function in_quarantine(locdat, p, leakage = 0.05)::Bool
    @assert 0.0 <= leakage <= 1.0 "leakage value must be between 0.0 and 1.0 inclusive"
    if leakage == 1.0
        return false  # everyone leaks quarantine is always false
    elseif leakage == 0.0
        return locdat.quar[p] # no one leaks--depends on the person's status
    else
        rand() < leakage ? false : locdat.quar[p]
    end
end