####################################################
# progression.jl for ilm model
#     how folks who have gotten sick progress through stages of the disease in the simulation:
#           progression
####################################################

    
"""
    progression!(p, infectset, progressionset, vaxset, dovax, riskshift!, probvec, locdat)

People who have become infected progress through conditions from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. Finally,  
they move to recovered or dead.
"""
@inline function progression!(locdat, p, infectset, progressionset, vaxset, dovax, probvec)

    today = DAY_CTR[:day]

    @inbounds person = locdat[p]

    # extract traits for a single person p where p is the row index in locdat
    p_duration = person.duration  # no. of days p has been infected
    p_cond = person.cond         # p's condition
    p_status = person.status     # p's status, etc.
    p_agegrp = person.agegrp  
    p_vaxstatus = person.vaxstatus
    p_recovday = p_status === :recovered ? person.recovday[end] : 0
    p_variant = person.variant[end]

    probtree = progressionset[p_variant].tree   

    # if person's agegrp and duration match a progression stage, 
    #     get the progression array for rows=from condition and cols=to condition/status
    pr_arr = get(getfield(probtree, Symbol(p_agegrp)),  # getfield returns dict for an agegrp
                p_duration,                             # key of the agegrp dict: duration of disease at progression checkpoint
                [])                                     # if criteria not met: default to empty

    if isempty(pr_arr)  
        person.duration += 1  # no progression today: one more day with current condition
    else
        # person p may progress to another condition of the disease, recover or die
        probvec[:] = pr_arr[mapcondition(p_cond), :] # probabilities of recovery, nil, mild, sick, severe, dead given p's current condition

        # effect on severity and progressing based on recovery from a previous infection
        recoveff =  if p_status == :recovered
                        recoveffect(today, p_recovday, p_variant, infectset)
                    else
                        1.0
                    end
                    
        # effect on severity and progressing based on being vaccinated
        vaxeff = if p_vaxstatus === :none
                     1.0
                 else
                     p_vaxrcvd = person.vaxrcvd[end]
                     p_vaxday = person.vaxday[end]
                     p_variant = person.variant[end]
                     vaxeffect(today, infectset, vaxset, p_vaxstatus, p_variant, p_vaxrcvd, 
                                p_vaxday, mode=:progression)
                 end

        risk = riskfactor(recoveff, vaxeff)
        
        if dovax     # vaccination changes probability, thus timing, of progressing to different stages of disease
            redistribute_probability!(probvec, risk, p_duration) 
        end

        doprogression!(person, probvec)   
    end
    
    return nothing 
end


function riskfactor(recoveff, vaxeff)
    combinedfactor = min(recoveff, vaxeff)
    riskfactor = clamp(combinedfactor, 0.0, 0.97)  
end


"""
Redistribute the probabilty of progressing through conditions based on
vaccination, recovery from prior infection and the variant of the patient.
Modifies in place the probvec: probalities of progressing to each state
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
    doprogression!(progression_cols, p, probvec)

Progress an infected person to a new condition or status if called
with a progression vector (trvec) or increment
the number of days the person has been sick.
"""
function doprogression!(person, probvec)

    outcome = categorical_sim(probvec) # which outcome based on probability...?

    tocond = map_progression(outcome) # 1->recover, 2->nil, 3->mild, 4->sick, 5->severe, 6->dead  see data_mapping.jl

    if tocond == :dead  
        person.deadday = DAY_CTR[:day]
        person.status = :dead  # change the status
        person.cond = :uninfected # change the condition   TODO: could we leave this as it was to summarize how people died
    elseif tocond == :recovered
        push!(person.recovday, DAY_CTR[:day])
        person.status = :recovered
        person.cond = :uninfected   # TODO decide if this makes sense--using this to maintain a history of past infection
    else  # tocond is another disease condition 
        person.cond = tocond   # change the condition = degree of sickness
        person.duration += 1    # advance number of days person has been sick
    end    
end
