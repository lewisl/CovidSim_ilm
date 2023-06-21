####################################################
# progression.jl for ilm model
#     how folks who have gotten sick progress through stages of the disease in the simulation:
#           progression
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
@inline function progression!(p, infectset, progressionset, vaxset, dovax, probvec, prog_cols)

    today = DAY_CTR[:day]

    # extract traits for a single person p where p is the row index in locdat
    @inbounds begin
        p_duration = prog_cols.duration[p]  # no. of days p has been infected
        p_cond = prog_cols.cond[p]          # p's condition
        p_status = prog_cols.status[p]      # p's status, etc.
        p_agegrp = prog_cols.agegrp[p]  
        p_vaxstatus = prog_cols.vaxstatus[p]
        p_recovday = prog_cols.recovday[p][end]
        p_variant = prog_cols.variant[p][end]
    end

    prtree = progressionset[p_variant].tree   

    # if person's agegrp and duration match a progression stage, get the progression array for rows=from and cols=to
    pr_arr = get( 
                getfield(prtree, Symbol(p_agegrp)),  # tree for an agegrp
                p_duration,                          # duration of disease: day on which condition may progress
                [])                                  # if selection criteria not met:  empty

    if !isempty(pr_arr)  # person p may progress to another condition of the disease, recover, or die
        probvec[:] = pr_arr[mapcondition(p_cond), :] # probabilities of recovery, nil, mild, sick, severe, dead given p's current condition

        # effect on severity and progressioning based on recovery from a previous infection
        recoveff =  @inbounds if p_status == :recovered
                        recoveffect(today, p_recovday, p_variant, infectset)
                    else
                        1.0
                    end

        vaxeff = @inbounds if p_vaxstatus === :none
                        1.0
                    else
                        p_vaxrcvd = prog_cols.vaxrcvd[p][end]
                        p_vaxday = prog_cols.vaxday[p][end]
                        p_variant = prog_cols.variant[p][end]
                        # effect on severity and progressing based on being vaccinated
                        vaxeffect(today, infectset, vaxset, p_vaxstatus, p_variant, p_vaxrcvd, p_vaxday, mode=:progression)
                    end

        risk = riskfactor(recoveff, vaxeff)
        
        if dovax
            redistribute_probability!(probvec, risk, p_duration) 
        end

        doprogression!(prog_cols, p, probvec)

    else # no progression today!
        c_duration[p] += 1  # one more day in current condition
    end
    
    return    
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
function doprogression!(prog_cols, p, probvec)

    choice = categorical_sim(probvec) # which outcome based on probability...?

    # debugging
    @assert choice != 0 "choice of to condition resulted in 0. Must be 1 through 6"

    tocond = map_progression(choice) # 1->recover, 2->nil, 3->mild, 4->sick, 5->severe, 6->dead  see data_mapping.jl

    if tocond == :dead  
        @inbounds begin
        prog_cols.deadday[p] = DAY_CTR[:day]
        prog_cols.status[p] = :dead  # change the status
        c_cond[p] = :uninfected # change the condition
        end
    elseif tocond == :recovered
        @inbounds begin
        push!(prog_cols.recovday[p], DAY_CTR[:day])
        prog_cols.status[p] = :recovered
        progression_cols.cond[p] = :uninfected   # TODO decide if this makes sense--using this to maintain a history of past infection
        end
    else  # tocond to another infectious condition 
        @inbounds begin
        prog_cols.cond[p] = tocond   # change the condition = degree of sickness
        prog_cols.duration[p] += 1    # advance number of days person has been sick
        end
    end    
end
