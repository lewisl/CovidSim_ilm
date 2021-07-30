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

# %% [markdown]
# ## Quick Test of Group Model ###

# %%
using Markdown
using InteractiveUtils

# %%
cd(joinpath(homedir(),"Dropbox/Online Coursework/Covid/group-src"))

# %%
using CovidSim_group

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :default, agegrps)

# %%
result_dict, series = run_a_sim(180, 38015, showr0=false, silent=true, runcases=[seed_1_6]);

# %%
cumplot(series, 38015)

# %%
keys(result_dict)

# %%
result_dict["dat"]["openmx"][38015]

# %%
