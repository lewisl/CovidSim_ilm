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
@inline function transition!(locdat, infect_idx, transitionset, dovax, dovariant)
        

    # aliases for person attribute columns--deref the named tuple once
    v_sickday = locdat.sickday
    v_cond = locdat.cond
    v_agegrp = locdat.agegrp

    if dovax == true
    elseif dovariant == true
    else
        dectree = transitionset[:new]
    end

    for p in infect_idx  # p for infected person    
        p_sickday = v_sickday[p]
        p_cond = v_cond[p]
        p_agegrp = v_agegrp[p]  # agegroup of person p = agegrp column of locale data, row p 

        # if person's agegrp, sickday, and condition match a transition stage
        transvec = has(dectree[p_agegrp], p_sickday, p_cond) 

        dotransition!(locdat, p, p_cond, transvec) # perform transition logic and update population table

    end  
end


function has(agetr::Dict, sickday::Int, p_cond::condition)::Union{Vector{Float64}, Nothing}
    for (stage, v) in agetr
        if v[:sickday] == sickday
            trvec = collect(v[:transition][p_cond, :])
            if sum(trvec) > 0.0
                return trvec
            end
        end
    end
    return nothing
end


"""
    dotransition!(locdat, p, p_cond, trvec::Union{Vector{Float64}, Nothing})

Transition an infected person to a new condition or status if called
with a transition vector (trvec) or increment
the number of days the person has been sick.
"""
@inline function dotransition!(locdat, p, p_cond, trvec::Union{Vector{Float64}, Nothing})
   
    if isnothing(trvec)
        # if locdat.sickday[p] >= 25
        #     println("$(day_ctr[:day]): agegrp: $(locdat.agegrp[p]) sickday: $(locdat.sickday[p]) from cond: $p_cond")
        # end

        locdat.sickday[p] += 1  

    else
        choice = shift(categorical_sim(trvec)) # which outcome...?

        tocond = conditionshift(p_cond, choice)  # next condition or status

        # if locdat.sickday[p] >= 25
        #     println("$(day_ctr[:day]): agegrp: $(locdat.agegrp[p]) sickday: $(locdat.sickday[p]) from cond: $p_cond to cond: $tocond")
        # end

        if tocond == dead  
            locdat.deadday[p] = day_ctr[:day]
            locdat.status[p] = dead  # change the status
            locdat.cond[p] = notsick # change the condition
        elseif tocond == recovered
            locdat.recovday[p] = day_ctr[:day]
            locdat.status[p] = recovered
            locdat.cond[p] = notsick
        else   
            locdat.cond[p] = tocond   # change the condition = degree of sickness
            locdat.sickday[p] += 1    # advance number of days person has been sick
        end    


    end

end


function conditionshift(now, sh)

    if sh == same
            now
    elseif sh == die
        dead
    elseif sh == recover
        recovered
    elseif sh == improve
        Int(now) - 1 < Int(typemin(condition)) ? typemin(condition) : condition(Int(now) - 1)
    elseif sh == worse
        Int(now) + 1 > Int(typemax(condition)) ? typemax(condition) : condition(Int(now) + 1)
    elseif sh == worseplus
        Int(now) + 2 > Int(typemax(condition)) ? typemax(condition) : condition(Int(now) + 2)
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
function travelin!(dat=popdat)
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