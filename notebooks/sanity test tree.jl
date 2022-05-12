# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Julia 1.7.2
#     language: julia
#     name: julia-1.7
# ---

# %%
using CovidSim_ilm

# %%
cs = CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML
using PrettyPrint
using LinearAlgebra
using Dates

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# # Test model setup

# %% tags=[]
ndays = 180
locale = 38015
model180 = buildsim(ndays, locale;  
    day1 = Date("2020-01-01", "yyyy-mm-dd"),
    dovax = true,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantfilename = "variants.yml",
);

# %%
keys(model180)


# ### Social Parameters

# %%
fieldnames(typeof(model180.social))

# %%
typeof(model180.social.gammashape)

# %%
model180.social.contactfactors 

# %%
model180.social.contactfactors[Int(sick)-4, Int(age0_19)]

# %%
touchfactors = model180.social.touchfactors

# %% [markdown]
# ## Infectset

# %%
model180.infectset

# %% [markdown]
# ### Parameters for infection based on variant

# %%
fieldnames(typeof(model180.infectset[:omicron_ba1]))

# %%
println(CovidSim_ilm.variantlist)

# %%
pprintln(keys(model180.infectset))

# %%
pprintln(model180.infectset[:base])

# %%
sendbase = model180.infectset[:base].sendrisk
recvbase = model180.infectset[:base].recvrisk

# %%
pprint(model180.infectset[:alpha])


# %% [markdown]
# ## Transition parameters

# %%
model180.transitionset # the decision transition matrices for all age groups are loaded

# %%
basetransition = model180.transitionset[:base]
println("Types...")
println(typeof(model180.transitionset))
println(typeof(basetransition))
println(typeof(basetransition.tree))
println(typeof(basetransition.tree.age0_19))
println(typeof(basetransition.tree.age0_19[5]))
println()
pprintln(basetransition)

# %%
fieldnames(typeof(basetransition.tree))

# %%
keys(basetransition.tree.age0_19)

# %%
age0_19tree = basetransition.tree.age0_19

# %%
age0_19tree[5]

# %%
basetransition.factors

# %%
basetransition.factors.vaxhalflifeadjust


# %% [markdown]
# ## Design a new sanity check

# %% [markdown]
# ### Check rows

# %%
anywrong = false
for checkpt in [5,9,14,19,25]
    for row in 1:4
        if isapprox(sum(age0_19tree[checkpt][row,:]), 1.0)
        else
            println("row probabilities not equal 1.0: $checkpoint, $row: $(age0_19[checkpt][row,:])")
            anywrong = true
        end
    end
end
anywrong || println("All probability rows sum to 1.0")

# %% [markdown]
# ### Check outcomes are all dead or recovered

# %%
transitionset = model180.transitionset
variant = :base
setage = age0_19
breakday = 5
checkpoints = [5,9,14,19,25]
agetrees = transitionset[variant].tree

fieldnames(typeof(agetrees))

# %%
thisage = age80_up
for checkday in [5,9,14,19,25]
    println("******* checkday $checkday")
    display(getfield(agetrees, Symbol(thisage))[checkday])
    println()
end

# %%
torecovered = Dict(k => 0.0 for k in [nil, mild, sick, severe])
todead = Dict(k => 0.0 for k in [nil, mild, sick, severe])

# %%
currage = age0_19
thisagetree = getfield(agetrees, Symbol(currage))
res = copy(thisagetree[5])
firstpass = true
for checkpt in checkpoints[2:end]
    println("*** $checkpt")
    # torecovered[nil] += sum(view(res, :, 1))
    # todead[nil] += sum(view(res, :, 6))
    # display(res)
    

    trconds = transpose(copy(res[:, 2:5]))
    if firstpass
        for i in 1
            res[:] = trconds[:, i] .* thisagetree[checkpt]
            torecovered[nil] += sum(res[:, 1])
            todead[nil] += sum(res[:, 6])
            # trconds = transpose(copy(res[:, 2:5]))
        end
        firstpass = false
    else
        for i in 1:4
            res[:] = trconds[:, i] .* thisagetree[checkpt]
            torecovered[nil] += sum(res[:, 1])
            todead[nil] += sum(res[:, 6])
            # trconds = transpose(copy(res[:, 2:5]))
        end
    end

    println("recovered ", torecovered[nil], " dead ", todead[nil], " all ", torecovered[nil]+todead[nil])
end


# %%
thisagetree[25]

# %%
res

# %%
trconds

# %%
age0_19tree[9]

# %%
trconds[:,1] .* age0_19tree[9]

# %%
trconds[:,2] .* age0_19tree[9]

# %%
currage = age0_19
thisagetree = getfield(agetrees, Symbol(currage))

# %%
cs.sanitycheck(agetrees)

# %%
