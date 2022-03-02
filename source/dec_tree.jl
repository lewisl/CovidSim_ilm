
#############################################################
# dec_tree.jl
# decision tree for transition
#############################################################


        Base.@kwdef struct Transitiondef
            sickday::Int64
            transition::Matrix{Float64}
        end

Base.@kwdef struct Agetree
    age0_19::Vector{Transitiondef}
    age20_39::Vector{Transitiondef}
    age40_59::Vector{Transitiondef}
    age60_79::Vector{Transitiondef}
    age80_up::Vector{Transitiondef}
end


Base.@kwdef struct Transitionfactors
    riskadjust::Union{Vector{Float64}, Nothing}
    vaxhalflifeadjust::Union{Dict{Symbol, Float64}, Nothing}
    
        # inner method
        function Transitionfactors(factordict)
            riskadj = get(factordict, "riskadjust", nothing)
            vaxadj = get(factordict, "vaxhalflifeadjust", nothing)
            vaxadj = if !isnothing(vaxadj)
                        Dict(Symbol(k)=>v for (k,v) in vaxadj)
                     end
            new(riskadj, vaxadj)
        end
end

Base.@kwdef struct Transitionparams
    tree::Union{Agetree, Nothing}
    factors::Transitionfactors   # use [] for nothing
end


# method for creating from file per variant
function setup_dt(dtfilename::String)
    arrays = YAML.load_file(dtfilename)

    newdict = (
        Dict(symbol2agegrp(Symbol(k1)) =>      # agegrp    
            Dict(k2 =>                    # stage in 1:5
                Dict(Symbol(k3) =>  if k3 == :sickday
                                        v3
                                    else 
                                        Dict(symbol2condition(k4) => v4 for (k4, v4) in v3)
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
                    newdict[k1][k2][k3] = OffsetArray(out, 5:8, 1:6)   # index by condition from nil to severe, 1:6
                end
            end
        end
    end

    return newdict
end

# method for creating from Dict to nested structs
"""
    function setup_dt(trdict::Dict)

A decision tree is of the following type and structure:
- an instance of struct Agetree with members:
    - age0_19
    - age20_39
    - age40_59
    - age60_79
    - age80_up
- Each age field of Agetree is a vector of struct Transitiondef with fields:
    - sickday: the day on which transitions to different disease outcomes occurs
    - transition: an array that maps from current conditions (rows) to outcomes (columns that are the probability of
     either a different sickness condition or a final outcome of either recover or dead) 
another disease state (condition or status).

"""
function setup_dt(trdict::Dict)

    prepdict = Dict(age_key => [Transitiondef(brk[:sickday], 
                                    vcat(brk[:transition][:nil]',brk[:transition][:mild]',
                                        brk[:transition][:sick]',brk[:transition][:severe]')) 
                                for (_, brk) in sort(age)] 
                        for (age_key, age) in sort(trdict))

    return Agetree([prepdict[age] for age in keys(sort(trdict))]...) # sort and splat the arguments without using keywords
end


"""
transitionT is Type alias for type that holds a transition decision tree
"""
const transitionT = Dict{agegrp, Dict{Int64, Dict{Symbol, Any}}}



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

function display_tree_array(tree)
    for agegrp in keys(tree)
        agetree = tree[agegrp]
        println("agegrp: ", agegrp, " =>")
        for brkday_idx in keys(agetree)
            sickdaytree = agetree[brkday_idx]
            println("    sickday: ", sickdaytree[:sickday])
            println("    transitions: ")
            for r in eachrow(sickdaytree[:transition])
                print("      "); println(r)
            end
        end  # for sickday
    end   # for agegrp     
end


function display_tree_struct(tree)
    for agegrp in fieldnames(typeof(tree))
        agetree = getfield(tree, Symbol(agegrp))
        println(agegrp, " # field of struct Agetree, values are vectors")
        for brk in 1:length(agetree)
            println("    brk: $brk", " # element of Vector{Transitiondef}")
            println("    sickday: ", agetree[brk].sickday)
            println("    transitions: ")
            for r in eachrow(agetree[brk].transition)
                print("        "); println(r)
            end
        end
    end
end

function sanitycheck(dectree::Dict)
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


function sanitycheck(dectree::Agetree)
    for age in fieldnames(Agetree)
        seqs = getseqs(getfield(dectree, Symbol(age)))
        probs, allpr = verifyprobs(seqs)
        println("for agegroup ", age)
        for p in pairs(probs)
            println("    ",p)
        end
        println("    Prob total: ",allpr)
    end
end

"""
Use for Dict representation of trees. Find all sequences of conditions by transition date and current condition through to new conditions
for a single agegrp.
"""
function getseqs(dt_this_age::Dict; maxsearches = 100)
    # find the top nodes
    dt_this_age = sort(dt_this_age)
    breakdays = collect(keys(dt_this_age))
    k1 = first(breakdays)
    todo = [] # array of node sequences 
    done = [] # ditto

    # gather the outcomes at the first breakday for the starting conditions
    # no transition has happened yet: these are initial conditions: the first sequence(s) to be extended
    for fromcond in keys(dt_this_age[k1])
        for i in 1:length(dt_this_age[k1][fromcond][:outcomes])
            outcome = dt_this_age[k1][fromcond][:outcomes][i]
            prob = dt_this_age[k1][fromcond][:probs][i]
            push!(todo, [(sickday=k1, fromcond=fromcond, tocond=outcome, prob=prob)])
        end
    end


    # build sequences from top to terminal states: recovered or dead
    ctr = 0
    while !isempty(todo)
        if (ctr += 1) > maxsearches
            @assert false "maxsearches exceeded when buildig sequences through transition tree at $ctr"
        end
        seq = popfirst!(todo)  
        lastnode = seq[end]
        breakday, fromcond, tocond = lastnode
        nxtidx = findfirst(isequal(breakday), breakdays) + 1
        for brk in breakdays[nxtidx:end]
            if tocond in keys(dt_this_age[brk])   # keys are the fromcond at the next break day so previous tocond == current fromcond
                for i in 1:length(dt_this_age[brk][tocond][:outcomes])
                    outcome = dt_this_age[brk][tocond][:outcomes][i]
                    prob = dt_this_age[brk][tocond][:probs][i]
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



"""
Use for Array and struct representation of trees. Find all sequences of conditions by transition date and current condition 
through to new conditions for a single agegrp.
"""
function getseqs(dt_age::Vector{Transitiondef}; maxsearches = 100)
    # dt_age is an vector of Transitiondef
    # find the top nodes (node is a condition/status transition)
    # dt_age = sort(dt_age)  # in order by breaks
    breakdays = [brk.sickday for brk in dt_age]
    brk_idx = 1:size(dt_age,1)
    todo = [] # array of node sequences 
    done = [] # ditto

    # gather the outcomes at the first breakday for the starting conditions
    # no transition has happened yet: these are initial conditions: the first sequence(s) to be extended

    breakday = breakdays[1]                      # dt_age[brk1][:sickday]
    transitions = dt_age[1].transition  # from-to matrix of probabilities

    for row in eachindex(transitions[:,1])     
        for i in 1:size(transitions,2)              # to outcome -> column within row
            if transitions[row, i] != 0.0       
                outcome = transition_cases[i]
                prob = transitions[row, i] 
                push!(todo, [(sickday=breakday, fromcond=map2cond(row), tocond=outcome, prob=prob)])
            end
        end
    end

    # build sequences from top to terminal states: recovered or dead
    ctr = 0
    while !isempty(todo)
        if (ctr += 1) > maxsearches
            @assert false "maxsearches exceeded when buildig sequences through transition tree at $ctr"
        end
        seq = popfirst!(todo)  
        lastnode = seq[end]
        breakday, fromcond, tocond = lastnode    # breakday = day of transition; fromcond = row index; tocond = column index
        nxtidx = findfirst(isequal(breakday), breakdays) + 1
        for brk in brk_idx[nxtidx:end]
            next_transition = dt_age[brk].transition

            breakday = breakdays[brk]
            next_steps = [ i for i in eachindex(next_transition[:,1]) if any(next_transition[i,:] .!= 0.0) ] 

            if map2cond(tocond) in next_steps # keys are the fromcond at the next break day so previous tocond == current fromcond
                outcomes_idx = findall(next_transition[map2cond(tocond),:] .!= 0.0)  

                for i in 1:length(outcomes_idx)                      # 1:length(dt_age[brk][tocond][:outcomes])
                    outcome = outcomes_idx[i]  # dt_age[brk][tocond][:outcomes][i]
                    prob = next_transition[map2cond(tocond), outcome]    # prob = dt_age[brk][tocond][:probs][i]   
                    newseq = vcat(seq, (sickday=breakday, fromcond=tocond, tocond=transition_cases[outcome], prob=prob))
                    outcome = transition_cases[outcome]
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

# what a tree looks like using nested structs and from->to array for transitions
#=
age0_19 =               # field of struct Agetree, value is a vector of Transitiondef, sorted by sickday
    [   sickday: 5      # field of Transitiondef
        transitions:    # field of Transitiondef
            [0.0, 0.391304347826087, 0.4891304347826087, 0.11956521739130435, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 9
        transitions: 
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.95, 0.05, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 14
        transitions: 
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.85, 0.0, 0.0, 0.12, 0.03, 0.0]
            [0.692, 0.0, 0.0, 0.0, 0.302, 0.006]
        sickday: 19
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.891, 0.0, 0.0, 0.0, 0.106, 0.003]
        sickday: 25
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.976, 0.0, 0.0, 0.0, 0.0, 0.024]
            [0.91, 0.0, 0.0, 0.0, 0.0, 0.09]
    ]
age20_39 =
    [   sickday: 5
        transitions: 
            [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 9
        transitions: 
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            [0.85, 0.0, 0.0, 0.15, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 14
        transitions: 
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.83, 0.0, 0.0, 0.1, 0.07, 0.0]
            [0.474, 0.0, 0.0, 0.0, 0.514, 0.012]
        sickday: 19
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
        sickday: 25
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
            [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
    ]
age40_59 =
    [   sickday: 5
        transitions: 
            [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 9
        transitions: 
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            [0.85, 0.0, 0.05, 0.1, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 14
        transitions: 
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            [0.85, 0.0, 0.0, 0.14, 0.01, 0.0]
            [0.776, 0.0, 0.0, 0.0, 0.206, 0.018]
        sickday: 19
        transitions: 
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
        sickday: 25
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
            [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
    ]
age60_79 =
        transitions: 
            [0.0, 0.15, 0.6, 0.25, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 9
        transitions: 
            [0.62, 0.0, 0.0, 0.38, 0.0, 0.0]
            [0.5, 0.0, 0.25, 0.25, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.78, 0.22, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 14
        transitions: 
            [0.8, 0.1, 0.1, 0.0, 0.0, 0.0]
            [0.8, 0.0, 0.15, 0.05, 0.0, 0.0]
            [0.8, 0.0, 0.0, 0.1, 0.1, 0.0]
            [0.165, 0.0, 0.0, 0.0, 0.715, 0.12]
        sickday: 19
        transitions: 
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
        sickday: 25
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.76, 0.0, 0.0, 0.0, 0.0, 0.24]
            [0.688, 0.0, 0.0, 0.0, 0.0, 0.312]
    ]
age80_up =
    [   sickday: 5
        transitions: 
            [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 9
        transitions: 
            [0.5, 0.0, 0.0, 0.5, 0.0, 0.0]
            [0.0, 0.0, 0.4, 0.6, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.6, 0.4, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        sickday: 14
        transitions: 
            [0.7, 0.0, 0.3, 0.0, 0.0, 0.0]
            [0.7, 0.0, 0.0, 0.3, 0.0, 0.0]
            [0.7, 0.0, 0.0, 0.1, 0.2, 0.0]
            [0.12, 0.0, 0.0, 0.0, 0.67, 0.21]
        sickday: 19
        transitions: 
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.49, 0.0, 0.0, 0.0, 0.24, 0.27]
        sickday: 25
        transitions: 
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.682, 0.0, 0.0, 0.0, 0.0, 0.318]
            [0.676, 0.0, 0.0, 0.0, 0.0, 0.324]
    ]
=#



# what an older version tree built as a Dict looks like using from->to array for transitions
#=
agegrp: age0_19 =>
    sickday: 25
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.976, 0.0, 0.0, 0.0, 0.0, 0.024]
      [0.91, 0.0, 0.0, 0.0, 0.0, 0.09]
    sickday: 19
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.891, 0.0, 0.0, 0.0, 0.106, 0.003]
    sickday: 9
    transitions: 
      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.95, 0.05, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sickday: 14
    transitions: 
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.85, 0.0, 0.0, 0.12, 0.03, 0.0]
      [0.692, 0.0, 0.0, 0.0, 0.302, 0.006]
    sickday: 5
    transitions: 
      [0.0, 0.4, 0.5, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age40_59 =>
    sickday: 25
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
      [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
    sickday: 19
    transitions: 
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
    sickday: 9
    transitions: 
      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.85, 0.0, 0.05, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sickday: 14
    transitions: 
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.85, 0.0, 0.0, 0.14, 0.01, 0.0]
      [0.776, 0.0, 0.0, 0.0, 0.206, 0.018]
    sickday: 5
    transitions: 
      [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age20_39 =>
    sickday: 25
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
      [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
    sickday: 19
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
    sickday: 9
    transitions: 
      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.85, 0.0, 0.0, 0.15, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sickday: 14
    transitions: 
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.83, 0.0, 0.0, 0.1, 0.07, 0.0]
      [0.474, 0.0, 0.0, 0.0, 0.514, 0.012]
    sickday: 5
    transitions: 
      [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age60_79 =>
    sickday: 25
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.76, 0.0, 0.0, 0.0, 0.0, 0.24]
      [0.688, 0.0, 0.0, 0.0, 0.0, 0.312]
    sickday: 19
    transitions: 
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
    sickday: 9
    transitions: 
      [0.62, 0.0, 0.0, 0.38, 0.0, 0.0]
      [0.5, 0.0, 0.25, 0.25, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.78, 0.22, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sickday: 14
    transitions: 
      [0.8, 0.1, 0.1, 0.0, 0.0, 0.0]
      [0.8, 0.0, 0.15, 0.05, 0.0, 0.0]
      [0.8, 0.0, 0.0, 0.1, 0.1, 0.0]
      [0.165, 0.0, 0.0, 0.0, 0.715, 0.12]
    sickday: 5
    transitions: 
      [0.0, 0.15, 0.6, 0.25, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age80_up =>
    sickday: 25
    transitions: 
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.682, 0.0, 0.0, 0.0, 0.0, 0.318]
      [0.676, 0.0, 0.0, 0.0, 0.0, 0.324]
    sickday: 19
    transitions: 
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.49, 0.0, 0.0, 0.0, 0.24, 0.27]
    sickday: 9
    transitions: 
      [0.5, 0.0, 0.0, 0.5, 0.0, 0.0]
      [0.0, 0.0, 0.4, 0.6, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.6, 0.4, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sickday: 14
    transitions: 
      [0.7, 0.0, 0.3, 0.0, 0.0, 0.0]
      [0.7, 0.0, 0.0, 0.3, 0.0, 0.0]
      [0.7, 0.0, 0.0, 0.1, 0.2, 0.0]
      [0.12, 0.0, 0.0, 0.0, 0.67, 0.21]
    sickday: 5
    transitions: 
      [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

=#




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


