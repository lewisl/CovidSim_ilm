####################################################
# transition.jl for ilm model
#     change status of folks in simulation:
#           transition
#           travel
####################################################

Base.@kwdef struct Transitiondef
    sickday::Int64
    trarr::Matrix{Float64}
end

Base.@kwdef struct Agetree
    age0_19::Vector{Transitiondef}
    age20_39::Vector{Transitiondef}
    age40_59::Vector{Transitiondef}
    age60_79::Vector{Transitiondef}
    age80_up::Vector{Transitiondef}
end

Base.@kwdef struct Transitionparams
    tree::Union{Agetree, Nothing}
    riskadjust::Vector{Float64}   # use [] for nothing
end


    
"""
    transition!(locdat, infect_idx, dectree)

People who have become infectious transition through cases from
nil (asymptomatic) to mild to sick to severe, depending on their
agegroup, days of being exposed, and some probability. Finally,  
they move to recovered or dead.

locdat must be a population table for a single locale.
"""
@inline function transition!(locdat, infect_idx, spreadset, transitionset, vaxset, dovax, dovariant; vaxfn! = noop, trvec = zeros(6))
        
    # aliases for person attribute columns--deref the named tuple once
    v_sickday = locdat.sickday
    v_cond = locdat.cond
    v_agegrp = locdat.agegrp

    if dovax
        vaxfn! = vaxfn! == noop ? vaxtransitioneffect! : vaxfn!  # last branch new vaxfn! was passed in
    end


    transarray = transitionset[:base]   # TODO test variant of each person

    for p in infect_idx  # p for infected person    
        p_sickday = v_sickday[p]
        p_cond = v_cond[p]
        p_agegrp = v_agegrp[p]  # agegroup of person p = agegrp column of locale data, row p 

        # if person's agegrp, sickday, and condition match a transition stage
        transvec = has(transarray[p_agegrp], p_sickday, p_cond, trvec) 

        riskadjustments = riskadjust(spreadset, vaxset, locdat, p)

        vaxfn!(transvec, riskadjustments) #vaxfn! will be function noop or function vaxtransitioneffect

        dotransition!(locdat, p, p_cond, transvec) # perform transition logic and update population table

    end  
end


function has(agetr::Dict, sickday::Int, p_cond::condition, trvec)::Union{Vector{Float64}, Nothing}
    for (stage, v) in agetr
        if v[:sickday] == sickday
            trvec[:] = v[:transition][map2cond(p_cond), :]
            if sum(trvec) > 0.0
                return trvec
            end
        end
    end
    return nothing
end


# this will be set to the variable vaxfn!
function vaxtransitioneffect!(transvec, adjustments)
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
            # locdat.cond[p] = uninfected
        else   
            locdat.cond[p] = tocond   # change the condition = degree of sickness
            locdat.sickday[p] += 1    # advance number of days person has been sick
        end    


    end

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