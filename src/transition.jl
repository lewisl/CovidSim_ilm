####################################################
# transition.jl for ilm model
#     change status of folks in simulation:
#           transition
#           travel
####################################################

    
"""
    transition!(p, infectset, transitionset, vaxset, dovax, vaxfn!, transvec, <columns of locdat>)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. Finally,  
they move to recovered or dead.

Required columns of locdat are cond, status, agegrp, duration, sdcomply, variant, vaxstatus, recovday
vaxrcvd, vaxday, deadday.
"""
@inline function transition!(p, infectset, transitionset, vaxset, dovax, vaxfn!, transvec,
            c_cond,
            c_status,
            c_agegrp,
            c_duration,
            c_sdcomply,
            c_variant,
            c_vaxstatus,
            c_recovday,
            c_vaxrcvd,
            c_vaxday,
            c_deadday
            )

    if dovax
        vaxfn! = vaxfn! == noop ? vaxtransitioneffect! : vaxfn!  # last branch new vaxfn! was passed in
    end

    # extract traits for this person p where p is the row index in locdat
    @inbounds begin
        p_duration = c_duration[p]
        p_cond = c_cond[p]
        p_status = c_status[p]
        p_agegrp = c_agegrp[p]  
        p_vaxstatus = c_vaxstatus[p]
        p_recovday = c_recovday[p][end]
        p_variant = c_variant[p][end]
    end

    trtree = transitionset[p_variant].tree   

    # if person's agegrp and duration match a transition stage
    tr_arr = get( getfield(trtree, Symbol(p_agegrp)), p_duration, [])

    if !isempty(tr_arr)  # let's transition person p 
        transvec[:] = tr_arr[mapcondition(p_cond), :] # probabilities of recovery, nil, mild, sick, severe, dead given current condition

        # effect on severity and transitioning based on recovery from a previous infection
        recoveff =  @inbounds if p_status == recovered
                        tr_recoveffect(p_recovday, p_variant, infectset)
                    else
                        1.0
                    end

        # effect on severity and transitioning based on being vaccinated
        vaxeff = @inbounds if p_vaxstatus != :none
                        p_vaxrcvd = c_vaxrcvd[p][end]
                        p_vaxday = c_vaxday[p][end]
                        p_variant = c_variant[p][end]
                        tr_vaxeffect(infectset, vaxset, p_vaxstatus, p_variant, p_vaxrcvd, p_vaxday)
                    else
                        1.0
                    end
        
        vaxfn!(transvec, recoveff, vaxeff) #vaxfn! will be function noop or function vaxtransitioneffect

        dotransition!(p, transvec, # perform transition logic and update population table->must pass columns, not scalars  
                        c_duration,
                        c_deadday,
                        c_status,
                        c_cond,
                        c_recovday
                    )
    else # no transition
        c_duration[p] += 1  # one more day in current condition
    end
    
    return    
end


# this function set to the variable vaxfn!
@inline function vaxtransitioneffect!(transvec, recoveff, vaxeff, varianteff=1.0)
    @inbounds begin
    combinedfactor = varianteff * min(recoveff, vaxeff)
    riskfactor = squashfunc(combinedfactor)  

    for c in (sick, severe, dead)
        transvec[maptransition(c)] *= riskfactor
    end

    correction = 1.0 / sum(transvec)  # normalize to sum to 1.0
    transvec[:] .*= correction
    end
end

"""
    dotransition!(locdat, p, p_cond, trvec::Union{Vector{Float64}, Nothing})

Transition an infected person to a new condition or status if called
with a transition vector (trvec) or increment
the number of days the person has been sick.
"""
function dotransition!(p, transvec,         
                c_duration,
                c_deadday,
                c_status,
                c_cond,
                c_recovday
            )

        choice = categorical_sim(transvec) # which outcome based on probability...?

        # debugging
        @assert choice != 0 "choice of to condition resulted in 0. Must be 1 through 6"

        tocond = maptransition(choice) # 1->recover, 2->nil, 3->mild, 4->sick, 5->severe, 6->dead  see data_mapping.jl

        if tocond == dead  
            @inbounds begin
            c_deadday[p] = day_ctr[:day]
            c_status[p] = dead  # change the status
            c_cond[p] = uninfected # change the condition
            end
        elseif tocond == recovered
            @inbounds begin
            push!(c_recovday[p], day_ctr[:day])
            c_status[p] = recovered
            c_cond[p] = uninfected   # TODO decide if this makes sense--using this to maintain a history of past infection
            end
        else  # tocond to another infectious condition 
            @inbounds begin
            c_cond[p] = tocond   # change the condition = degree of sickness
            c_duration[p] += 1    # advance number of days person has been sick
            end
        end    
    # end
end


"""
    vaximmunity(infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday)

Immunity from vaccination for a single person.
"""
@inline function tr_vaxeffect(infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday)

    today = day_ctr[:day]
    oneshotfactor = 0.85   # TODO yet another parameter to put somewhere...!

    # person's vaccine conditions
    days_after_vax = today - vaxday
    @assert today >= days_after_vax "today's date must be >= to day of most recent shot"

    # vaccine characteristics
    halflife = vaxset[vaxrcvd].halflife
    full_effect_days = vaxset[vaxrcvd].full_effect_days
    infectfactor = vaxset[vaxrcvd].infectfactor
    tr_vaxeffect = vaxset[vaxrcvd].infectreduce[vaxstatus][spr_variant]

    # rise and decay of vaccine effectiveness
    rise_lower=0.5    # TODO need to make these inputs somewhere
    decay_lower=0.10   # lindecay argument
    csig = 10.0       # sigdecay argument
    rise = riseup(today - days_after_vax, full_effect_days, rise_lower, 1.0)
    days_after_full_effect = today - (days_after_vax + full_effect_days)     #clamp(today - (lastshotday + full_effect_days), 0, Int)
    # decay = lindecay(days_after_full_effect, halflife, decay_lower)
    decay = sigdecay(days_after_full_effect, halflife, csig=csig)
    vaxmod = rise * decay

    factor = 1.0 - (vaxmod * tr_vaxeffect * infectfactor) 

    return factor
end


"""
    tr_recoveffect(recovday, targ_variant, spr_variant, infectset)

Immunity from recovery for a single person.
"""
@inline function tr_recoveffect(recovday, targ_variant, infectset)

        today = day_ctr[:day]
        days_post_recov = today - recovday 

        if days_post_recov > 0   # TODO should be an assert: does this run day of or day after recovery?
            # get the max immunity
            immstrength = infectset[targ_variant].recovery_immunity[targ_variant]

            # get the declined value
            decay_lower=0.1   # lindecay argument
            csig = 10.0       # sigdecay argument
            immhalflife = infectset[targ_variant].immunehalflife
            # immdecline = lindecay(days_post_recov, immhalflife, decay_lower)
            immdecline = sigdecay(days_post_recov, immhalflife, csig=csig)

            factor = 1.0 - (immdecline * immstrength)
        else
            factor = 1.0
        end

    return factor
end


"""
For a locale, randomly choose the number of people from each agegroup with
condition of {unexposed, infectious, recovered} who travel to each
other locale. Add to the travelq.
"""
function travelout!(fromloc, locales, rules=[])    # TODO THIS WON'T WORK ANY MORE!
    # 10.5 microseconds for 5 locales
    # choose distribution of people traveling by age and condition:
        # unexposed, infectious, recovered -> ignore duration for now
    # TODO: more frequent travel to and from Major and Large cities
    # TODO: should the caller do the loop across locales?   YES
    travdests = collect(locales)
    deleteat!(travdests,findfirst(isequal(fromloc), travdests))
    bins = lim = length(travdests) + 1
    for agegrp in agegrps
        for cond in [unexposed, infectious, recovered]
            name = string(cond)
            for duration in durations
                numfolks = sum(grab(cond, agegrp, duration, fromloc)) # the from locale, all durations
                travcnt = floor(Int, gamma_prob(travprobs[agegrp]) * numfolks)  # interpret as fraction of people who will travel
                x = rand(travdests, travcnt)  # randomize across destinations
                bydest = bucket(x, vals=1:length(travdests))
                for dest in 1:length(bydest)
                    isempty(bydest) && continue
                    cnt = bydest[dest]
                    iszero(cnt) && continue
                    enqueue!(travelq, travitem(cnt, fromloc, dest, agegrp, duration, name))
                end
            end
        end
    end
end


"""
Assuming a daily cycle, at the beginning of the day
process the queue of travelers from the end of the previous day.
Remove groups of travelers by agegrp, duration, and condition
from where they departed.  Add them to their destination.
"""
function travelin!(dat=popdat)   # TODO THIS DOESN'T WORK ANYMORE
    while !isempty(travelq)
        g = dequeue!(travelq)
        cond = eval(Symbol(g.cond))
        minus!(g.cnt, cond, g.agegrp, g.duration, g.from, dat=dat)
        plus!(g.cnt, cond, g.agegrp, g.duration, g.to, dat=dat)
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