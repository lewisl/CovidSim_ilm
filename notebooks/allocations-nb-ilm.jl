# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %%
using CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
using PrettyPrint

# %%
cd(joinpath(homedir(),"Dropbox/Online Coursework/Covid/ilm-src"))

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :default, agegrps)

# %% [markdown]
# # Run a simulation for 40 days

# %%
result_dict, series = run_a_sim(40, 38015, showr0=false, silent=true, runcases=[seed_1_6]);

# %%
keys(result_dict)

# %%
dectrees = result_dict["dt_dict"]

# %%
popdat = result_dict["dat"]["popdat"]


# %% [markdown]
# ### Check allocations of transition! function

# %%
@time transition!(popdat, 38015, dectrees)

# %% [markdown]
# (2.46 k allocations: 87.141 KiB): based on number of times through the loop
#
# what's the cost per iteration?

# %%
dectree = dectrees["dt"] # before starting the loop--once per function call

# %%
@time dectree_a1 = dectree[1]  # within the loop; no allocations ~ 5 ns

# %% jupyter={"outputs_hidden": true} tags=[]
pprint(dectree_a1)

# %%
node = dectree[1][14][7] # once per iteration

# %%
@btime $dectree[1][14][7]["probs"]

# %% [markdown]
# ##### allocations: (1 allocation: 32 bytes) or 0 allocs: 0 bytes if empty

# %%
typeof(node)

# %%
@time node["probs"]

# %%
@time node["outcomes"]

# %% [markdown]
# #### transition! is not bad for allocations, but there are a lot of iterations in a simulation
# It would be a fair amount of work to eliminate the allocation completely.

# %% [markdown]
# ## try a different data structure

# %%
holder = Dict(14 => Dict(6=>Dict(:probs=>[1.0], :outcomes=>[3]),
                         7=>Dict(:probs=>[0.85, 0.12, 0.03], :outcomes=>[3,7,8]),
                         8=>Dict(:probs=>[0.692, 0.302, 0.006], :outcomes=>[3,8,4])), 
               9 => Dict(5=>Dict(:probs=>[.9, .1], :outcomes=>[3,7]),
                         6=>Dict(:probs=>[1.0], :outcomes=>[6]),
                         7=>Dict(:probs=>[.95, .05], :outcomes=>[7,8])))

# %%
choice = CovidSim_ilm.categorical_sim(holder[14][7][:probs]) # no allocations here

# %%
tocond = holder[14][7][:outcomes][choice] # no allocation here

# %%
yaml_str = "1:                                         
  5:
    5:                                     
      probs: [0.4, 0.5, 0.1]
      outcomes: [5, 6, 7]
  9:
    5:                                     
      probs: [0.9, 0.1]
      outcomes: [0,7]
    6:   
      probs: [1.0]
      outcomes: [6]                                 
    7:   
      probs: [0.95, 0.05]                                   
      outcomes: [7,8]"

# %%
using YAML

# %%
YAML.load(yaml_str)

# %% [markdown]
# ## Create trees from YAML file with different structure

# %%
dtfilename="../parameters/dec_tree_all_25.yml"

# %%
trees = YAML.load_file(dtfilename)
    # next: change 2nd level keys from 2 item array{Int} [9, 5] to Tuple{Int, Int} (9,5)
    trees = Dict(i => Dict(Tuple(k)=>trees[i][k] for k in keys(trees[i])) for i in keys(trees))

    # next: change the type of next node item from array{Int} [25, 8] to Tuple{Int, Int} (25, 8)
    for agegrp in agegrps
        for (k,v) in trees[agegrp]
           for item in v
                item["next"] = Tuple(item["next"])
            end
        end
    end

# %% jupyter={"outputs_hidden": true} tags=[]
pprint(trees)

# %%
newdict = Dict()
for agegrp in agegrps
    newdict[agegrp] = Dict()
    for node in keys(trees[agegrp])
        a = node[1]
        b = node[2]
            probs = [branch["pr"] for branch in trees[agegrp][node]]
            outcomes = [branch["tocond"] for branch in trees[agegrp][node]]
            branches = [branch for branch in trees[agegrp][node]]
        if haskey(newdict[agegrp], a)
            newdict[agegrp][a][b] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
        else
            newdict[agegrp][a]=Dict()
            newdict[agegrp][a][b] = Dict("probs"=>probs, "outcomes"=>outcomes, "branches"=>branches)
        end
    end
end
newdict = Dict(i=>sort(newdict[i], rev=true) for i in agegrps) 

# %%
newdict[5]

# %%
newdict[5][14]

# %%
newdict[5][14][6]

# %%
@btime $newdict[5][14][6]["probs"]

# %% [markdown]
# test it...

# %%
categorical_sim = CovidSim_ilm.categorical_sim

# %%
function partial_transition(tree, agegrp, sickday_p, cond_p, foo=0)
    # node = get_node2(tree[agegrp], sickday_p, cond_p)
    # if isempty(node) # no transition for this person based on sickday and condition
    if !haskey(tree, sickday_p) || !haskey(tree[sickday_p], cond_p)
        # @assert p_tup.sickday < sickdaylim "Person made it to last day and was not removed:\n     $p_tup\n"
        foo += 1
    else  # change of the person p's state--a transition
        node = tree[agegrp][sickday_p][cond_p]
        choice = categorical_sim(node["probs"]) # rand(Categorical(node["probs"])) # which branch...?
        tocond = node["outcomes"][choice]
        return tocond
    end
end

# %%
agegrp = 1; sickday_p = 14; cond_p = 6

# %%
@btime partial_transition($newdict, $agegrp, $sickday_p, $cond_p)

# %%
180 * 1000 * (9 * 1e-9)

# %%
isa(6, Dict)

# %%
@btime newdict[agegrp][sickday_p][cond_p];

# %% [markdown]
# ### Check allocations of spreading (2 functions)

# %%
@time spread!(popdat, 38015, [], spreadparams) # returns (n_spreaders, n_contacts, n_touched, n_newly_infected)

# %% [markdown]
# not terrible: 0.005223 seconds (1.29 k allocations: 899.457 KiB, 69.25% compilation time)
#
# we could pre-allocate spread_idx and contactable_idx at the pop size. would need signal value for the effective
# end of the valid part of the array.

# %%
