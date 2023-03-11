####################################################
# progression.jl for ilm model
#     change condition or status of folks who have gotten sick in the simulation:
#           progression
#     travel
####################################################

    
"""
    progression!(p, infectset, progressionset, vaxset, dovax, riskshift!, probvec, <columns of locdat>)

People who have become infectious progress through conditions from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. Finally,  
they move to recovered or dead.

Required columns of locdat are cond, status, agegrp, duration, sdcase, variant, vaxstatus, recovday,
vaxrcvd, vaxday, deadday.
"""
@inline function progression!(p, infectset, progressionset, vaxset, dovax, probvec,
            c_cond,
            c_status,
            c_agegrp,
            c_duration,
            c_sdcase,
            c_variant,
            c_vaxstatus,
            c_recovday,
            c_vaxrcvd,
            c_vaxday,
            c_deadday
            )

    today = DAY_CTR[:day]


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

    prtree = progressionset[p_variant].tree   

    # if person's agegrp and duration match a progression stage, get the progression array for rows=from and cols=to
    pr_arr = get( getfield(prtree, Symbol(p_agegrp)), p_duration, [])

    if !isempty(pr_arr)  # let's progress person p 
        probvec[:] = pr_arr[mapcondition(p_cond), :] # probabilities of recovery, nil, mild, sick, severe, dead given current condition

        # effect on severity and progressioning based on recovery from a previous infection
        recoveff =  @inbounds if p_status == :recovered
                        recoveffect(today, p_recovday, p_variant, infectset)
                    else
                        1.0
                    end

        # effect on severity and progressioning based on being vaccinated
        vaxeff = @inbounds if p_vaxstatus === :none
                        1.0
                    else
                        p_vaxrcvd = c_vaxrcvd[p][end]
                        p_vaxday = c_vaxday[p][end]
                        p_variant = c_variant[p][end]
                        vaxeffect(today, infectset, vaxset, p_vaxstatus, p_variant, p_vaxrcvd, p_vaxday, mode=:progression)
                    end

        risk = riskfactor(recoveff, vaxeff)
        
        if dovax
            redistribute_probability!(probvec, risk, p_duration) 
        end

        doprogression!(p, probvec, # perform progression logic and update population table->must pass columns, not scalars  
                        c_duration,
                        c_deadday,
                        c_status,
                        c_cond,
                        c_recovday
                    )
    else # no progression
        c_duration[p] += 1  # one more day in current condition
    end
    
    return    
end


function riskfactor(recoveff, vaxeff)
    combinedfactor = min(recoveff, vaxeff)
    riskfactor = squashfunc(combinedfactor)  
end


"""
Redistribute the probabilty of progressing through conditions based on
vaccination, recovery from prior infection and the variant of the patient.
"""
@inline function redistribute_probability!(probvec, riskfactor, duration)
    @inbounds begin

        excessprob = 0.0
        for toprob in (:sick, :severe, :dead)  
            idx = map_progression(toprob) # in file data_mapping.jl: map Symbol of condition or status to integer index in the probability vector
            excess1 = probvec[idx] * (1.0 - riskfactor)
            probvec[idx] = probvec[idx] - excess1  # reduce probability of serious outcomes
            excessprob += excess1
        end

        if duration == DURATIONLIM   # clear anyone left to recovered or dead
            idx = map_progression(:recovered)
            probvec[idx] = probvec[idx] + excessprob
        else
            excessprob = excessprob / 3.0
            for toprob in (:recovered, :nil, :mild)
                idx = map_progression(toprob)
                probvec[idx] = probvec[idx] + excessprob  # redistribute probability to less serious outcomes
            end
        end

    end
end


"""
    doprogression!(locdat, p, p_cond, trvec::Union{Vector{Float64}, Nothing})

Progress an infected person to a new condition or status if called
with a progression vector (trvec) or increment
the number of days the person has been sick.
"""
function doprogression!(p, probvec,         
                c_duration,
                c_deadday,
                c_status,
                c_cond,
                c_recovday
            )

        choice = categorical_sim(probvec) # which outcome based on probability...?

        # debugging
        @assert choice != 0 "choice of to condition resulted in 0. Must be 1 through 6"

        tocond = map_progression(choice) # 1->recover, 2->nil, 3->mild, 4->sick, 5->severe, 6->dead  see data_mapping.jl

        if tocond == :dead  
            @inbounds begin
            c_deadday[p] = DAY_CTR[:day]
            c_status[p] = :dead  # change the status
            c_cond[p] = :uninfected # change the condition
            end
        elseif tocond == :recovered
            @inbounds begin
            push!(c_recovday[p], DAY_CTR[:day])
            c_status[p] = :recovered
            c_cond[p] = :uninfected   # TODO decide if this makes sense--using this to maintain a history of past infection
            end
        else  # tocond to another infectious condition 
            @inbounds begin
            c_cond[p] = tocond   # change the condition = degree of sickness
            c_duration[p] += 1    # advance number of days person has been sick
            end
        end    
    # end
end



# TODO   This is ancient code and WILL NOT WORK in current version of simulaton code

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
    for agegrp in AGEGRPS
        for cond in [unexposed, infectious, recovered]
            name = string(cond)
            for duration in DURATIONS
                numfolks = sum(grab(cond, agegrp, duration, fromloc)) # the from locale, all DURATIONS
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