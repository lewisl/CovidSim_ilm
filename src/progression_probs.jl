
######################################################################
# progression_probs.jl
#   Create probabilities of progressing through stages of the illness,
#   based on age group and variant of the virus.
#   Estimates are adjusted to result in death rates roughly in line 
#   with published death rates by age group.
#
# Set up decision matrix based on input file
# Provide utilities for checking plausibility of progression probabilities
######################################################################


"""
Hold progression matrices by age group.
"""
Base.@kwdef struct Agetree
  age0_19::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
  age20_39::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
  age40_59::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
  age60_79::Dict{Int, Matrix{Float64}} =  Dict{Int, Matrix{Float64}}()
  age80_up::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
end

"""
Hold factors to alter progression matrices off base for
different vaccines and variants.
"""
Base.@kwdef struct ProgressionFactors
  riskadjust::Union{Vector{Float64}, Nothing}
  vaxhalflifeadjust::Union{Dict{Symbol, Float64}, Nothing}
  
      # inner method
      function ProgressionFactors(factordict)
          riskadj = get(factordict, :riskadjust, nothing)
          vaxadj = get(factordict, :vaxhalflifeadjust, nothing)
          vaxadj = if !isnothing(vaxadj)
                      Dict(Symbol(k)=>v for (k,v) in vaxadj)
                   end
          new(riskadj, vaxadj)
      end
end

Base.@kwdef struct ProgressionParams
  tree::Union{Agetree, Nothing}
  factors::ProgressionFactors   # use [] for nothing
end


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
    prepdict = Dict(age_key => Dict(parse(Int, string(duration)) =>
                                    vcat(prob_vecs[:nil]',  prob_vecs[:mild]',  
                                         prob_vecs[:sick]', prob_vecs[:severe]')
                                for (duration, prob_vecs) in sort(duration))  
                        for (age_key, duration) in sort(trdict))
    
    # top-level key of dict, agegrp, becomes a field value of the struct Agetree
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


function display_tree(tree::Agetree)
    conds = [:nil, :mild, :sick, :severe]
    for agegrp in fieldnames(typeof(tree))
        agedict = getfield(tree, Symbol(agegrp))
        println(agegrp, " # field of struct Agetree, which is Dict{Int, Matrix{Float64}}")
        for duration in sort(collect(keys(agedict)))
            println("  duration: ", duration, " = number of days being sick")
            println("    progression probabilities: ")
            println("    from    to")
            println("            recovered nil mild sick severe dead")
            for (i, r) in enumerate(eachrow(agedict[duration]))
                @printf("    %-7s ", conds[i]); println(r)
            end
            println()
        end
    end
end


function sanitycheck(tree::Agetree)
    for age in fieldnames(Agetree)
        println("\nfor agegroup ", age); flush(stdout)
        seqs = getseqs(getfield(tree, Symbol(age)))
        probs, allpr, restable = verifyprobs(seqs)
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

    # gather the outcomes at the first breakday for the starting conditions
    # no progression has happened yet: these are initial conditions: the first sequence(s) to be extended
    for fromcond in mapcondition(:nil)  # everyone starts at nil
        for i in 1:size(dt_this_age[k1],2)  # no. of columns
            outcome = map_progression(i)
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
                outcome = map_progression(i)
                prob = dt_this_age[brk][mapcondition(tocond), map_progression(outcome)]
                # @show(outcome, prob); println();
                newseq = vcat(seq, (duration=brk, fromcond=tocond, tocond=outcome, prob=prob))
                # @show(newseq); println()
                if prob > 0.0
                    if (outcome == :dead) | (outcome == :recovered)  # terminal node reached--no more nodes to add
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
    ret = Dict(:dead=>0.0, :recovered=>0.0)
    restable = Table(duration=[5,9,14,19,25], from=[:nil, :nil, :nil, :nil, :nil], recovered=[0.0,0.0,0.0,0.0,0.0], dead=[0.0,0.0,0.0,0.0,0.0])

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
    push!(restable, (duration=100, from=:uninfected, recovered=sum(restable.recovered), dead=sum(restable.dead)))
    return ret, allpr, restable
end

# This is the YAML input for a transtion tree, which is provided for the :base variant and optionally other variants.
# If a variant doesn't provide its own progression tree, it can adjust the :base tree.
#=
      age0_19:       
        1:
          duration: 5
          progression:
            nil:       [0.0, 0.4, 0.5, 0.1, 0.0, 0.0]
            mild:      [0.0, 0.3, 0.6, 0.1, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        2:
          duration: 9
          progression:
            nil:       [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            mild:      [0.4, 0.0, 0.6, 0.0, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 0.95, 0.05, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        3:
          duration: 14
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.85, 0.0, 0.0, 0.12, 0.03, 0.0]
            severe:    [0.692, 0.0, 0.0, 0.0, 0.302, 0.006]
        4:
          duration: 19
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            severe:    [0.891, 0.0, 0.0, 0.0, 0.106, 0.003]
        5:
          duration: 25
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.976, 0.0, 0.0, 0.0, 0.0, 0.024]
            severe:    [0.91, 0.0, 0.0, 0.0, 0.0, 0.09]
      age20_39:
        1:
          duration: 5
          progression:
            nil:       [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
            mild:      [0.0, 0.15, 0.85, 0.0, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            severe:    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        2:
          duration: 9
          progression:
            nil:       [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            mild:      [.85, 0.0, 0.0, 0.15, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        3:
          duration: 14
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.83, 0.0, 0.0, 0.1, 0.07, 0.0]
            severe:    [0.474, 0.0, 0.0, 0.0, 0.514, 0.012]
        4:
          duration: 19
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
            severe:    [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
        5:
          duration: 25
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
            severe:    [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
      age40_59:
        1:
          duration: 5
          progression:
            nil:       [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
            mild:      [0.0, 0.15, 0.85, 0.0, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        2:
          duration: 9
          progression:
            nil:       [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            mild:      [0.85, 0.0, .05, 0.1, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        3:
          duration: 14
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
            sick:      [0.85, 0.0, 0.0, 0.14, 0.01,  0.0]
            severe:    [0.776, 0.0, 0.0, 0.0, 0.206, 0.018]
        4:
          duration: 19
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
            severe:    [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
        5:
          duration: 25
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
            severe:    [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
      age60_79:
        1:
          duration: 5
          progression:
            nil:       [0.0, 0.15, 0.6, 0.25, 0.0, 0.0]
            mild:      [0.0, 0.0, 0.7, 0.3, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        2:
          duration: 9
          progression:
            nil:       [0.62, 0.0, 0.0, 0.38, 0.0, 0.0]
            mild:      [0.5, 0.0, .25, 0.25, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 0.78, 0.22, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        3:
          duration: 14
          progression:
            nil:       [0.8, 0.1, 0.1, 0.0, 0.0, 0.0]
            mild:      [0.8, 0.0, 0.15, 0.05, 0.0, 0.0]
            sick:      [0.8, 0.0, 0.0, 0.1, 0.1, 0.0]
            severe:    [0.165, 0.0, 0.0, 0.0, 0.715, 0.12]
        4:
          duration: 19
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
            severe:    [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
        5:
          duration: 25
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.76, 0.0, 0.0, 0.0, 0.0, 0.24]
            severe:    [0.688, 0.0, 0.0, 0.0, 0.0, 0.312]
      age80_up:
        1:
          duration: 5
          progression:
            nil:       [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
            mild:      [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
            sick:      [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.4, 0.6, 0.0]
        2:
          duration: 9
          progression:
            nil:       [0.5, 0.0, 0.0, 0.5, 0.0, 0.0]
            mild:      [0.0, 0.0, 0.4, 0.6, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 0.6, 0.4, 0.0]
            severe:    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        3:
          duration: 14
          progression:
            nil:       [0.7, 0.0, 0.3, 0.0, 0.0, 0.0]
            mild:      [0.7, 0.0, 0.0, 0.3, 0.0, 0.0]
            sick:      [0.7, 0.0, 0.0, 0.1, 0.2, 0.0]
            severe:    [0.12, 0.0, 0.0, 0.0, 0.67, 0.21]
        4:
          duration: 19
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
            severe:    [0.49, 0.0, 0.0, 0.0, 0.24, 0.27]
        5:
          duration: 25
          progression:
            nil:       [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            mild:      [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            sick:      [0.682, 0.0, 0.0, 0.0, 0.0, 0.318]
            severe:    [0.676, 0.0, 0.0, 0.0, 0.0, 0.324]
=#

# This is converted to an internal Julia datastructure for use during the simulation.
# The tree is stored in dict progressionset[:base], which is a struct ProgressionParams. 
# The field tree of ProgressionParams is a struct Agetree. Each field of Agetree is an agegrp.
# The agegrp field's value is a Dict{Int, Matrix{Float64}}. The key Int value is the number of 
# days someone is sick (duration) when the person progresses to a different stage of the disease.
# The value matrix is 4x6.  Each row is the probability of progressing from a disease condition to
# recover, a condition (one of: nil, mild, sick, severe) or to dieing (dead). Each row
# of the matrix must sum to 1.0 to be valid probabilities.
# A formatted printout of the data structure is:

#=
age0_19 # field of struct Agetree, which is Dict{Int, Matrix{Float64}}
  duration: 5 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.0, 0.4, 0.5, 0.1, 0.0, 0.0]
    mild    [0.0, 0.3, 0.6, 0.1, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 9 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
    mild    [0.4, 0.0, 0.6, 0.0, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 0.95, 0.05, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 14 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.85, 0.0, 0.0, 0.12, 0.03, 0.0]
    severe  [0.692, 0.0, 0.0, 0.0, 0.302, 0.006]

  duration: 19 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    severe  [0.891, 0.0, 0.0, 0.0, 0.106, 0.003]

  duration: 25 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.976, 0.0, 0.0, 0.0, 0.0, 0.024]
    severe  [0.91, 0.0, 0.0, 0.0, 0.0, 0.09]

age20_39 # field of struct Agetree, which is Dict{Int, Matrix{Float64}}
  duration: 5 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
    mild    [0.0, 0.15, 0.85, 0.0, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    severe  [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]

  duration: 9 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
    mild    [0.85, 0.0, 0.0, 0.15, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 14 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.83, 0.0, 0.0, 0.1, 0.07, 0.0]
    severe  [0.474, 0.0, 0.0, 0.0, 0.514, 0.012]

  duration: 19 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]
    severe  [0.922, 0.0, 0.0, 0.0, 0.072, 0.006]

  duration: 25 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]
    severe  [0.964, 0.0, 0.0, 0.0, 0.0, 0.036]

age40_59 # field of struct Agetree, which is Dict{Int, Matrix{Float64}}
  duration: 5 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.0, 0.2, 0.7, 0.1, 0.0, 0.0]
    mild    [0.0, 0.15, 0.85, 0.0, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 9 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
    mild    [0.85, 0.0, 0.05, 0.1, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 0.9, 0.1, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 14 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [0.9, 0.0, 0.0, 0.1, 0.0, 0.0]
    sick    [0.85, 0.0, 0.0, 0.14, 0.01, 0.0]
    severe  [0.776, 0.0, 0.0, 0.0, 0.206, 0.018]

  duration: 19 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]
    severe  [0.856, 0.0, 0.0, 0.0, 0.126, 0.018]

  duration: 25 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]
    severe  [0.958, 0.0, 0.0, 0.0, 0.0, 0.042]

age60_79 # field of struct Agetree, which is Dict{Int, Matrix{Float64}}
  duration: 5 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.0, 0.15, 0.6, 0.25, 0.0, 0.0]
    mild    [0.0, 0.0, 0.7, 0.3, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 9 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.62, 0.0, 0.0, 0.38, 0.0, 0.0]
    mild    [0.5, 0.0, 0.25, 0.25, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 0.78, 0.22, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 14 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.8, 0.1, 0.1, 0.0, 0.0, 0.0]
    mild    [0.8, 0.0, 0.15, 0.05, 0.0, 0.0]
    sick    [0.8, 0.0, 0.0, 0.1, 0.1, 0.0]
    severe  [0.165, 0.0, 0.0, 0.0, 0.715, 0.12]

  duration: 19 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]
    severe  [0.81, 0.0, 0.0, 0.0, 0.13, 0.06]

  duration: 25 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.76, 0.0, 0.0, 0.0, 0.0, 0.24]
    severe  [0.688, 0.0, 0.0, 0.0, 0.0, 0.312]

age80_up # field of struct Agetree, which is Dict{Int, Matrix{Float64}}
  duration: 5 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
    mild    [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
    sick    [0.0, 0.1, 0.5, 0.4, 0.0, 0.0]
    severe  [0.0, 0.0, 0.0, 0.4, 0.6, 0.0]

  duration: 9 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.5, 0.0, 0.0, 0.5, 0.0, 0.0]
    mild    [0.0, 0.0, 0.4, 0.6, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 0.6, 0.4, 0.0]
    severe  [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]

  duration: 14 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [0.7, 0.0, 0.3, 0.0, 0.0, 0.0]
    mild    [0.7, 0.0, 0.0, 0.3, 0.0, 0.0]
    sick    [0.7, 0.0, 0.0, 0.1, 0.2, 0.0]
    severe  [0.12, 0.0, 0.0, 0.0, 0.67, 0.21]

  duration: 19 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    severe  [0.49, 0.0, 0.0, 0.0, 0.24, 0.27]

  duration: 25 = number of days being sick
    progression probabilities: 
    from    to
            recovered nil mild sick severe dead
    nil     [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    mild    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sick    [0.682, 0.0, 0.0, 0.0, 0.0, 0.318]
    severe  [0.676, 0.0, 0.0, 0.0, 0.0, 0.324]
=#
