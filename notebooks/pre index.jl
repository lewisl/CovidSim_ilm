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
using TypedTables
using StatsBase
using SplitApplyCombine
using BenchmarkTools

# %%
using CovidSim_ilm

# %%
cd(joinpath(homedir(), "Dropbox/Online Coursework/Covid/ilm-src"))

# %%
tst = CovidSim_ilm.pop_data(100_000, age_dist=age_dist, cols="track")

# %% jupyter={"outputs_hidden": true} tags=[]
tst.status[rand(1:100_000,5000)] .= infectious

# %%
countmap(tst.status)

# %%
@time in_idx = findall(tst.status .== 2);

# %%
@time in_view = mapview(==(2),tst.status);

# %%
@time tst.status[in_view];

# %%
@time tst.status[in_idx];

# %%
@time groupinds(tst.status)

# %%
sick_ages_idx = groupinds(tst.agegrp[tst.status .== 2])

# %%
maximum(sick_ages_idx[4])

# %%
age3_sick_idx = findall((tst.agegrp .== 3) .& (tst.status .== 2))

# %%
@btime $tst.status[52520], $tst.agegrp[52520]

# %%
@btime sickdict = Dict(i=>findall(($tst.agegrp .== i) .& ($tst.status .== 2)) for i in 1:5);

# %%
@btime $sickdict[1][50]

# %%
foo = @elapsed for i in eachindex(tst.status)
    if tst.status[i] == 2
        tst.duration[i] = 1
    end
end

# %%
bar = @elapsed begin
    infect_idx = CovidSim_ilm.optfindall(==(2), tst.status, 0.5)
    for p in infect_idx
        tst.duration[p] = 2
    end
end

# %% [markdown]
# ### sample vs rand

# %%
@btime sample($age3_sick_idx, 5, replace=true)

# %%
n = length(age3_sick_idx)
@btime begin
            who = rand(1:n,5)
            age3_sick_idx[who]
        end

# %% [markdown]
# ### We like sample

# %%
