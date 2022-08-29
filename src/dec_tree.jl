
#############################################################
# dec_tree.jl
# decision tree for progression
#############################################################



# method for creating from Dict to nested structs
"""
    function setup_dt(trdict::Dict)
A decision tree is:
- an instance of struct Agetree with members:
    - age0_19
    - age20_39
    - age40_59
    - age60_79
    - age80_up
- Each field value of Agetree is a dict:
    - key is duration: the day on which progression to different disease outcomes occur;
    - value is progression array: maps from current conditions (rows) to outcomes (columns) based on the probability of
      progressioning from the current condition to a different infectious condition or a final outcome of recover or dead.


    function setup_dt(basetree::Agetree, adjust::Vector{Float64})
This method builds the decision tree based on the progression array for the variant called :base (required input), using an adjustment
vector for another variant.

"""
function setup_dt(trdict::Dict)

    prepdict = Dict(age_key => Dict(brk[:duration] =>
                                    vcat(brk[:progression][:nil]',  brk[:progression][:mild]',
                                         brk[:progression][:sick]', brk[:progression][:severe]')
                                for (_, brk) in sort(params))
                        for (age_key, params) in sort(trdict))

    return Agetree([prepdict[age] for age in keys(sort(trdict))]...) # splats to fields of Agetree: a dict per agegrp

end


function setup_dt(basetree::Agetree, adjust::Vector{Float64})

    prepdict = Dict{Symbol, Dict{Int, Matrix{Float64}}}()

    for age in Symbol.(AGEGRPS)
        prepdict[age] = getfield(basetree, age)
        for arr in values(prepdict[age])  # duration==key::Int, arr==progression array 4 x 6
            for r in eachrow(arr)
                if sum(r) != 0.0
                    r[:] = r .* adjust
                    sum_r = sum(r)
                    r[:] = r ./ sum_r
                end
            end
        end
    end

    return Agetree([prepdict[age] for age in keys(sort(prepdict))]...)
end


function display_tree(tree)
    for agegrp in keys(tree)
        agetree = tree[agegrp]
        println("agegrp: ", agegrp, " =>")
        for duration in keys(agetree)
            durationtree = agetree[duration]
            println("    duration: ", duration, " =>")
            for fromcond in keys(durationtree)
                condtree = durationtree[fromcond]
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
        end  # for duration
    end   # for agegrp     
end

function display_tree_array(tree)
    for agegrp in keys(tree)
        agetree = tree[agegrp]
        println("agegrp: ", agegrp, " =>")
        for brkday_idx in keys(agetree)
            durationtree = agetree[brkday_idx]
            println("    duration: ", durationtree[:duration])
            println("    progression: ")
            for r in eachrow(durationtree[:progression])
                print("      "); println(r)
            end
        end  # for duration
    end   # for agegrp     
end


function display_tree_struct(tree)
    for agegrp in fieldnames(typeof(tree))
        agetree = getfield(tree, Symbol(agegrp))
        println(agegrp, " # field of struct Agetree, values are vectors")
        for brk in eachindex(agetree)
            println("    brk: $brk", " # sickness duration day on which progression occurs")
            println("    duration: ", agetree[brk].duration)
            println("    progressions: ")
            for r in eachrow(agetree[brk].progression)
                print("        "); println(r)
            end
        end
    end
end


function sanitycheck(dectree::Agetree)
    for age in fieldnames(Agetree)
        println("\nfor agegroup ", age); flush(stdout)
        seqs = getseqs(getfield(dectree, Symbol(age)))
        probs, allpr, restable = verifyprobs(seqs)
        # for p in pairs(probs)
        #     println("        ",p); 
        # end; println(); # flush(stdout); 
        # flag = isapprox(allpr, 1.0) ? "OK " : "BAD"
        # print("    $flag "); print("Prob total: ", allpr); println(); flush(stdout)
        display(restable)
        println("Recovered + Dead probability =  ", restable.dead[6] + restable.recovered[6])
    end
end

"""
Use for Dict representation of trees. Find all sequences of conditions by progression date and current condition through to new conditions
for a single agegrp.
"""
function getseqs(dt_this_age::Dict; maxsearches = 100)
    # find the top nodes
    dt_this_age = sort(dt_this_age)
    breakdays = collect(keys(dt_this_age))
    k1 = first(breakdays)
    todo = [] # array of node sequences 
    done = [] # ditto

    # println(dt_this_age)

    # gather the outcomes at the first breakday for the starting conditions
    # no progression has happened yet: these are initial conditions: the first sequence(s) to be extended
    for fromcond in mapcondition(nil)  # everyone starts at nhil
        for i in 1:size(dt_this_age[k1],2)  # no. of columns
            outcome = mapprogression(i)
            prob = dt_this_age[k1][fromcond, i]
            if (prob != 0.0)    
                push!(todo, [(duration=k1, fromcond=mapcondition(fromcond), tocond=outcome, prob=prob)])
            end
        end
    end

    # build sequences from top to terminal states: recovered or dead
    ctr = 0
    while !isempty(todo)
        if (ctr += 1) > maxsearches
            @assert false "maxsearches exceeded when building sequences through progression tree at $ctr"
        end
        seq = popfirst!(todo)  
        lastnode = seq[end]
        breakday, fromcond, tocond = lastnode
        nxtidx = findfirst(isequal(breakday), breakdays) + 1
        for brk in breakdays[nxtidx:end]
            # @show(brk); println()
  
            for i in 1:size(dt_this_age[k1],2)  # no. of columns
                outcome = mapprogression(i)
                prob = dt_this_age[brk][mapcondition(tocond), mapprogression(outcome)]
                # @show(outcome, prob); println();
                newseq = vcat(seq, (duration=brk, fromcond=tocond, tocond=outcome, prob=prob))
                # @show(newseq); println()
                if prob > 0.0
                    if (outcome == dead) | (outcome == recovered)  # terminal node reached--no more nodes to add
                        push!(done, newseq)
                        # verbose == true && begin; println(done); println(); end
                    else  # not at a terminal outcome: still more nodes to add
                        push!(todo, newseq)
                    end
                end
            end
            break # we found the tocond as a matching fromcond
            
        end
    end

    return done
end



function verifyprobs(seqs)
    ret = Dict(dead=>0.0, recovered=>0.0)
    restable = Table(duration=[5,9,14,19,25], from=[nil, nil, nil, nil, nil], recovered=[0.0,0.0,0.0,0.0,0.0], dead=[0.0,0.0,0.0,0.0,0.0])
    allpr = 0.0

    for seq in seqs
        # verbose && println(seq)
        pr = mapreduce(x->getindex(x,:prob), *, seq)
            outcome = last(seq).tocond
            ret[outcome] += pr
            allpr += pr
            lastcond=last(seq).fromcond
            atduration = last(seq).duration
            rowidx = findfirst(restable.duration .== atduration)
            getproperty(restable, Symbol(outcome))[rowidx] += pr
            getproperty(restable, :from)[rowidx] = lastcond
    end
    push!(restable, (duration=100, from=uninfected, recovered=sum(restable.recovered), dead=sum(restable.dead)))
    return ret, allpr, restable
end

# A transtion tree is provided for the :base variant and optionally other variants.
# If a variant doesn't provide its own progression tree, it can adjust the :base tree.
# The tree is stored in dict progressionset[:base], which is a struct ProgressionParams in 
# the field tree.  The tree is a struct Agetree.
# At progressionset[:base].tree you find...
#=
age0_19 =               # field of struct Agetree.  The value of this field is a Dict{Int, Matrix{Float64}}
    {5 =>      
            [0.0, 0.391304347826087, 0.4891304347826087, 0.11956521739130435, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    9 =>
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.95, 0.05, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    14 =>
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.85, 0.0, 0.0, 0.12, 0.03, 0.0
            [0.692, 0.0, 0.0, 0.0, 0.302, 0.006]
    19 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.891, 0.0, 0.0, 0.0, 0.106, 0.003]
    25 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.976, 0.0, 0.0, 0.0, 0.0, 0.024
            [0.91, 0.0, 0.0, 0.0, 0.0, 0.09]
    }
age20_39 =
    {5 =>
            [0.0, 0.2, 0.7, 0.1, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    9 =>
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0
            [0.85, 0.0, 0.0, 0.15, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.9, 0.1, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    14 =>
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.83, 0.0, 0.0, 0.1, 0.07, 0.0
            [0.474, 0.0, 0.0, 0.0, 0.514, 0.012]
    19 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
    25 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.964, 0.0, 0.0, 0.0, 0.0, 0.036
            [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
    }
age40_59 =
    5 =>
            [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    9 =>
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0
            [0.85, 0.0, 0.05, 0.1, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.9, 0.1, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    14 =>
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.9, 0.0, 0.0, 0.1, 0.0, 0.0
            [0.85, 0.0, 0.0, 0.14, 0.01, 0.0
            [0.776, 0.0, 0.0, 0.0, 0.206, 0.018]
    19 =>
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
    25 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.958, 0.0, 0.0, 0.0, 0.0, 0.042
            [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
    }
age60_79 =
    {5 =>
            [0.0, 0.15, 0.6, 0.25, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    9 =>
            [0.62, 0.0, 0.0, 0.38, 0.0, 0.0
            [0.5, 0.0, 0.25, 0.25, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.78, 0.22, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    14 =>
            [0.8, 0.1, 0.1, 0.0, 0.0, 0.0
            [0.8, 0.0, 0.15, 0.05, 0.0, 0.0
            [0.8, 0.0, 0.0, 0.1, 0.1, 0.0
            [0.165, 0.0, 0.0, 0.0, 0.715, 0.12]
    19 =>
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
    25 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            [0.76, 0.0, 0.0, 0.0, 0.0, 0.24
            [0.688, 0.0, 0.0, 0.0, 0.0, 0.312]
    }
age80_up =
    {5 =>

            [0.0, 0.1, 0.5, 0.4, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    9 =>
            [0.5, 0.0, 0.0, 0.5, 0.0, 0.0
             0.0, 0.0, 0.4, 0.6, 0.0, 0.0
             0.0, 0.0, 0.0, 0.6, 0.4, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    14 =>
            [0.7, 0.0, 0.3, 0.0, 0.0, 0.0
             0.7, 0.0, 0.0, 0.3, 0.0, 0.0
             0.7, 0.0, 0.0, 0.1, 0.2, 0.0
             0.12, 0.0, 0.0, 0.0, 0.67, 0.21]
    19 =>
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0
             1.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.49, 0.0, 0.0, 0.0, 0.24, 0.27]
    25 =>
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0
             0.682, 0.0, 0.0, 0.0, 0.0, 0.318
             0.676, 0.0, 0.0, 0.0, 0.0, 0.324]
    }
=#



# what an older version tree built as a Dict looks like using from->to array for progression
#=
agegrp: age0_19 =>
    duration: 25

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.976, 0.0, 0.0, 0.0, 0.0, 0.024]
      [0.91, 0.0, 0.0, 0.0, 0.0, 0.09]
    duration: 19

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.891, 0.0, 0.0, 0.0, 0.106, 0.003]
    duration: 9

      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.95, 0.05, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    duration: 14

      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.85, 0.0, 0.0, 0.12, 0.03, 0.0]
      [0.692, 0.0, 0.0, 0.0, 0.302, 0.006]
    duration: 5

      [0.0, 0.4, 0.5, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age40_59 =>
    duration: 25

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
      [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
    duration: 19

      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
    duration: 9

      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.85, 0.0, 0.05, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    duration: 14

      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.85, 0.0, 0.0, 0.14, 0.01, 0.0]
      [0.776, 0.0, 0.0, 0.0, 0.206, 0.018]
    duration: 5

      [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age20_39 =>
    duration: 25

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
      [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
    duration: 19

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
    duration: 9

      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
      [0.85, 0.0, 0.0, 0.15, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    duration: 14

      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.83, 0.0, 0.0, 0.1, 0.07, 0.0]
      [0.474, 0.0, 0.0, 0.0, 0.514, 0.012]
    duration: 5

      [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age60_79 =>
    duration: 25

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.76, 0.0, 0.0, 0.0, 0.0, 0.24]
      [0.688, 0.0, 0.0, 0.0, 0.0, 0.312]
    duration: 19

      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
    duration: 9

      [0.62, 0.0, 0.0, 0.38, 0.0, 0.0]
      [0.5, 0.0, 0.25, 0.25, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.78, 0.22, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    duration: 14

      [0.8, 0.1, 0.1, 0.0, 0.0, 0.0]
      [0.8, 0.0, 0.15, 0.05, 0.0, 0.0]
      [0.8, 0.0, 0.0, 0.1, 0.1, 0.0]
      [0.165, 0.0, 0.0, 0.0, 0.715, 0.12]
    duration: 5

      [0.0, 0.15, 0.6, 0.25, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
agegrp: age80_up =>
    duration: 25

      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.682, 0.0, 0.0, 0.0, 0.0, 0.318]
      [0.676, 0.0, 0.0, 0.0, 0.0, 0.324]
    duration: 19

      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.49, 0.0, 0.0, 0.0, 0.24, 0.27]
    duration: 9

      [0.5, 0.0, 0.0, 0.5, 0.0, 0.0]
      [0.0, 0.0, 0.4, 0.6, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.6, 0.4, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    duration: 14

      [0.7, 0.0, 0.3, 0.0, 0.0, 0.0]
      [0.7, 0.0, 0.0, 0.3, 0.0, 0.0]
      [0.7, 0.0, 0.0, 0.1, 0.2, 0.0]
      [0.12, 0.0, 0.0, 0.0, 0.67, 0.21]
    duration: 5

      [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

=#




#  what a tree looks like for 5 agegrps
#= 
agegrp: age0_19 =>                                      #"agegrp:" is not in the dict  agegrp value is an enum agegrp
    duration: 5 =>                                       #"duration:" is not in the dict
        fromcond: nil =>                                #"fromcond:" is not in the dict fromcond value is an enum condition
            probs: => [0.4, 0.5, 0.1]
            outcomes: => Condition[nil, mild, sick]
    duration: 25 =>
        fromcond: severe =>
            probs: => [0.91, 0.09]
            outcomes: => Status[recovered, dead]
        fromcond: sick =>
            probs: => [0.976, 0.024]
            outcomes: => Status[recovered, dead]
    duration: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Condition[mild]
        fromcond: nil =>
            probs: => [0.9, 0.1]
            outcomes: => Enum{Int32}[recovered, sick]   # this array contains enums from condition and status
        fromcond: sick =>
            probs: => [0.95, 0.05]
            outcomes: => Condition[sick, severe]
    duration: 14 =>
        fromcond: severe =>
            probs: => [0.692, 0.302, 0.006]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Status[recovered]
        fromcond: sick =>
            probs: => [0.85, 0.12, 0.03]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    duration: 19 =>
        fromcond: severe =>
            probs: => [0.891, 0.106, 0.003]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age60_79 =>
    duration: 5 =>
        fromcond: nil =>
            probs: => [0.15, 0.6, 0.25]
            outcomes: => condition[nil, mild, sick]
    duration: 25 =>
        fromcond: severe =>
            probs: => [0.688, 0.312]
            outcomes: => Status[recovered, dead]
        fromcond: sick =>
            probs: => [0.76, 0.24]
            outcomes: => Status[recovered, dead]
    duration: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Condition[mild]
        fromcond: nil =>
            probs: => [0.62, 0.38]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.78, 0.22]
            outcomes: => Condition[sick, severe]
    duration: 14 =>
        fromcond: severe =>
            probs: => [0.165, 0.715, 0.12]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Status[recovered]
        fromcond: sick =>
            probs: => [0.8, 0.1, 0.1]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    duration: 19 =>
        fromcond: severe =>
            probs: => [0.81, 0.13, 0.06]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age80_up =>
    duration: 5 =>
        fromcond: nil =>
            probs: => [0.1, 0.5, 0.4]
            outcomes: => Condition[nil, mild, sick]
    duration: 25 =>
        fromcond: severe =>
            probs: => [0.676, 0.324]
            outcomes: => Status[recovered, dead]
        fromcond: sick =>
            probs: => [0.682, 0.318]
            outcomes: => Status[recovered, dead]
    duration: 9 =>
        fromcond: mild =>
            probs: => [0.4, 0.6]
            outcomes: => Condition[mild, sick]
        fromcond: nil =>
            probs: => [0.5, 0.5]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.6, 0.4]
            outcomes: => Condition[sick, severe]
    duration: 14 =>
        fromcond: severe =>
            probs: => [0.12, 0.67, 0.21]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [0.7, 0.3]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.7, 0.1, 0.2]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    duration: 19 =>
        fromcond: severe =>
            probs: => [0.49, 0.24, 0.27]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age20_39 =>
    duration: 5 =>
        fromcond: nil =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => Condition[nil, mild, sick]
    duration: 25 =>
        fromcond: severe =>
            probs: => [0.964, 0.036]
            outcomes: => Status[recovered, dead]
        fromcond: sick =>
            probs: => [0.964, 0.036]
            outcomes: => Status[recovered, dead]
    duration: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Condition[mild]
        fromcond: nil =>
            probs: => [0.85, 0.15]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.9, 0.1]
            outcomes: => Condition[sick, severe]
    duration: 14 =>
        fromcond: severe =>
            probs: => [0.474, 0.514, 0.012]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Status[recovered]
        fromcond: sick =>
            probs: => [0.83, 0.1, 0.07]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    duration: 19 =>
        fromcond: severe =>
            probs: => [0.922, 0.072, 0.006]
            outcomes: => Enum{Int32}[recovered, severe, dead]
agegrp: age40_59 =>
    duration: 5 =>
        fromcond: nil =>
            probs: => [0.2, 0.7, 0.1]
            outcomes: => Condition[nil, mild, sick]
    duration: 25 =>
        fromcond: severe =>
            probs: => [0.958, 0.042]
            outcomes: => Status[recovered, dead]
        fromcond: sick =>
            probs: => [0.958, 0.042]
            outcomes: => Status[recovered, dead]
    duration: 9 =>
        fromcond: mild =>
            probs: => [1.0]
            outcomes: => Condition[mild]
        fromcond: nil =>
            probs: => [0.9, 0.1]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.9, 0.1]
            outcomes: => Condition[sick, severe]
    duration: 14 =>
        fromcond: severe =>
            probs: => [0.776, 0.206, 0.018]
            outcomes: => Enum{Int32}[recovered, severe, dead]
        fromcond: mild =>
            probs: => [0.9, 0.1]
            outcomes: => Enum{Int32}[recovered, sick]
        fromcond: sick =>
            probs: => [0.85, 0.14, 0.01]
            outcomes: => Enum{Int32}[recovered, sick, severe]
    duration: 19 =>
        fromcond: severe =>
            probs: => [0.856, 0.126, 0.018]
            outcomes: => Enum{Int32}[recovered, severe, dead]
=#


