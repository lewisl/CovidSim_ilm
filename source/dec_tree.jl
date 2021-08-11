
#############################################################
# dec_tree.jl
# decision tree for transition
#############################################################


function setup_dt(dtfilename)
    arrays = YAML.load_file(dtfilename)

    newdict = (
        Dict(symtoage[Symbol(k1)] =>      # agegrp    
            Dict(k2 =>                    # stage in 1:5
                Dict(Symbol(k3) =>  if k3 == "sickday"
                                        v3
                                    else 
                                        Dict(symtocond[Symbol(k4)] => v4 for (k4, v4) in v3)
                                    end
                    for (k3, v3) in v2)
                for (k2, v2) in v1)
            for (k1, v1) in arrays)
    )

    # change :transition value to an array
    for (k1, v1) in newdict         # k1 is agegrp
        for (k2, v2) in v1          # k2 is stage in 1:5
            for (k3, v3) in v2      # k3 is :sickday or :transition
                if k3 == :transition
                    out = vcat(v3[nil]',v3[mild]',v3[sick]',v3[severe]') # stack the vectors
                    newdict[k1][k2][k3] = OffsetArray(out, 5:8, 1:6)   # index by condition from nil to severe, shift from recover to die
                end
            end
        end
    end

    return newdict
end


"""
transitionT is Type alias for type that holds a transition decision tree
"""
const transitionT = Dict{agegrp, Dict{Int64, Dict{Symbol, Any}}}




function old_setup_dt(dtfilename)
    trees = YAML.load_file(dtfilename)

    newdict = (
                Dict(symtoage[Symbol(k1)] =>         
                    Dict(k2 =>             
                        Dict(symtocond[Symbol(k3)] => 
                            Dict(Symbol(k4) => v4 for (k4, v4) in v3) 
                                                        for (k3, v3) in v2)
                                                            for (k2, v2) in v1)
                                                                for (k1,v1) in trees)
                )

    # Convert values in :outcomes to Enum condition or status
    for (k1,v1) in newdict
        for (k2,v2) in v1
            for (k3,v3) in v2
                for (k4, v4) in v3
                    if k4 == :outcomes
                        outs = [symtoallconds[Symbol(out)]  for out in v4]
                        newdict[k1][k2][k3][k4] = outs
                    end
                end
            end
        end
    end

    return newdict
end


function display_tree(tree)
    for agegrp in keys(tree)
        agetree = tree[agegrp]
        println("agegrp: ", agegrp, " =>")
        for sickday in keys(agetree)
            sickdaytree = agetree[sickday]
            println("    sickday: ", sickday, " =>")
            for fromcond in keys(sickdaytree)
                condtree = sickdaytree[fromcond]
                println("        fromcond: ", fromcond, " =>")
                print("            probs: => ")
                println(condtree[:probs])
                #
                print("            outcomes: => ")
                println(condtree[:outcomes])
                #
                # println("            branches: =>")
                # for branch in keys(condtree["branches"])
                #     println("                ", condtree["branches"][branch])   
                # end
            end  # for fromcond
        end  # for sickday
    end   # for agegrp     
end


function sanitycheck(dectree)
    for age in agegrps
        seqs = getseqs(dectree[age])
        probs, allpr = verifyprobs(seqs)
        println("for agegroup ", age)
        for p in pairs(probs)
            println("    ",p)
        end
        println("    Prob total: ",allpr)
    end
end


"""
Find all sequences of conditions by transition date and current condition through to new conditions
for a single agegrp.
"""
function getseqs(dt_by_age)
    # find the top nodes
    dt_by_age = sort(dt_by_age)
    breakdays = collect(keys(dt_by_age))
    k1 = first(breakdays)
    todo = [] # array of node sequences 
    done = [] # ditto

    # gather the outcomes at the first breakday for the starting conditions
    # no transition has happened yet: these are initial conditions: the first sequence(s) to be extended
    for fromcond in keys(dt_by_age[k1])
        for i in 1:length(dt_by_age[k1][fromcond][:outcomes])
            outcome = dt_by_age[k1][fromcond][:outcomes][i]
            prob = dt_by_age[k1][fromcond][:probs][i]
            push!(todo, [(sickday=k1, fromcond=fromcond, tocond=outcome, prob=prob)])
        end
    end


    # build sequences from top to terminal states: recovered or dead
    while !isempty(todo)
        seq = popfirst!(todo)  
        lastnode = seq[end]
        breakday, fromcond, tocond = lastnode
        nxtidx = findfirst(isequal(breakday), breakdays) + 1
        for brk in breakdays[nxtidx:end]
            if tocond in keys(dt_by_age[brk])   # keys are the fromcond at the next break day so previous tocond == current fromcond
                for i in 1:length(dt_by_age[brk][tocond][:outcomes])
                    outcome = dt_by_age[brk][tocond][:outcomes][i]
                    prob = dt_by_age[brk][tocond][:probs][i]
                    newseq = vcat(seq, (sickday=brk, fromcond=tocond, tocond=outcome, prob=prob))
                    if (outcome == dead) | (outcome == recovered)  # terminal node reached--no more nodes to add
                        push!(done, newseq)
                    else  # not at a terminal outcome: still more nodes to add
                        push!(todo, newseq)
                    end
                end
                break # we found the tocond as a matching fromcond
            end
        end
    end

    return done
end


function verifyprobs(seqs)
    ret = Dict(dead=>0.0, recovered=>0.0)
    allpr = 0.0

    for seq in seqs
        pr = mapreduce(x->getindex(x,:prob), *, seq)
        outcome = last(seq).tocond
        ret[outcome] += pr
        allpr += pr
    end
    return ret, allpr
end






#  what a tree looks like for 5 agegrps
#= 
agegrp: age0_19 =>                                      #"agegrp:" is not in the dict  agegrp value is an enum agegrp
    sickday: 5 =>                                       #"sickday:" is not in the dict
        fromcond: nil =>                                #"fromcond:" is not in the dict fromcond value is an enum condition
            probs: => [0.4, 0.5, 0.1]
            outcomes: => condition[nil, mild, sick]
    sickday: 25 =>
        fromcond: severe =>
            probs: => [0.91, 0.09]
            outcomes: => status[recovered, dead]
        fromcond: sick =>
            probs: => [0.976, 0.024]
            outcomes: => status[recovered, dead]
    sickday: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => condition[mild]
        fromcond: nil =>
            probs: => [0.9, 0.1]
            outcomes: => Enum{Int32}[recovered, sick]   # this array contains enums from condition and status
        fromcond: sick =>
            probs: => [0.95, 0.05]
            outcomes: => condition[sick, severe]
    sickday: 14 =>
        fromcond: severe =>
            probs: => [0.692, 0.302, 0.006]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => status[recovered]
        fromcond: sick =>
            probs: => [0.85, 0.12, 0.03]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    sickday: 19 =>
        fromcond: severe =>
            probs: => [0.891, 0.106, 0.003]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age60_79 =>
    sickday: 5 =>
        fromcond: nil =>
            probs: => [0.15, 0.6, 0.25]
            outcomes: => condition[nil, mild, sick]
    sickday: 25 =>
        fromcond: severe =>
            probs: => [0.688, 0.312]
            outcomes: => status[recovered, dead]
        fromcond: sick =>
            probs: => [0.76, 0.24]
            outcomes: => status[recovered, dead]
    sickday: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => condition[mild]
        fromcond: nil =>
            probs: => [0.62, 0.38]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.78, 0.22]
            outcomes: => condition[sick, severe]
    sickday: 14 =>
        fromcond: severe =>
            probs: => [0.165, 0.715, 0.12]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => status[recovered]
        fromcond: sick =>
            probs: => [0.8, 0.1, 0.1]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    sickday: 19 =>
        fromcond: severe =>
            probs: => [0.81, 0.13, 0.06]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age80_up =>
    sickday: 5 =>
        fromcond: nil =>
            probs: => [0.1, 0.5, 0.4]
            outcomes: => condition[nil, mild, sick]
    sickday: 25 =>
        fromcond: severe =>
            probs: => [0.676, 0.324]
            outcomes: => status[recovered, dead]
        fromcond: sick =>
            probs: => [0.682, 0.318]
            outcomes: => status[recovered, dead]
    sickday: 9 =>
        fromcond: mild =>
            probs: => [0.4, 0.6]
            outcomes: => condition[mild, sick]
        fromcond: nil =>
            probs: => [0.5, 0.5]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.6, 0.4]
            outcomes: => condition[sick, severe]
    sickday: 14 =>
        fromcond: severe =>
            probs: => [0.12, 0.67, 0.21]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [0.7, 0.3]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.7, 0.1, 0.2]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    sickday: 19 =>
        fromcond: severe =>
            probs: => [0.49, 0.24, 0.27]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age20_39 =>
    sickday: 5 =>
        fromcond: nil =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => condition[nil, mild, sick]
    sickday: 25 =>
        fromcond: severe =>
            probs: => [0.964, 0.036]
            outcomes: => status[recovered, dead]
        fromcond: sick =>
            probs: => [0.964, 0.036]
            outcomes: => status[recovered, dead]
    sickday: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => condition[mild]
        fromcond: nil =>
            probs: => [0.85, 0.15]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.9, 0.1]
            outcomes: => condition[sick, severe]
    sickday: 14 =>
        fromcond: severe =>
            probs: => [0.474, 0.514, 0.012]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => status[recovered]
        fromcond: sick =>
            probs: => [0.83, 0.1, 0.07]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    sickday: 19 =>
        fromcond: severe =>
            probs: => [0.922, 0.072, 0.006]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age40_59 =>
    sickday: 5 =>
        fromcond: nil =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => condition[nil, mild, sick]
    sickday: 25 =>
        fromcond: severe =>
            probs: => [0.958, 0.042]
            outcomes: => status[recovered, dead]
        fromcond: sick =>
            probs: => [0.958, 0.042]
            outcomes: => status[recovered, dead]
    sickday: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => condition[mild]
        fromcond: nil =>
            probs: => [0.9, 0.1]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.9, 0.1]
            outcomes: => condition[sick, severe]
    sickday: 14 =>
        fromcond: severe =>
            probs: => [0.776, 0.206, 0.018]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [0.9, 0.1]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.85, 0.14, 0.01]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    sickday: 19 =>
        fromcond: severe =>
            probs: => [0.856, 0.126, 0.018]
            outcomes: => Enum{Int32}[recovered, severe, dead]
=#


