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
using BenchmarkTools
using Distributions

# %% [markdown]
# This was a good start on a new approach to transition. A few details had to be changed, combined
# but it did work first time algorithmically.

# %% [markdown]
# # Grab a population matrix after 80 days

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :base, agegrps)

# %%
alldat, series = run_a_sim(80, 38015, showr0=false, silent=true, runcases=[seed_1_6]);

# %%
locdat = alldat.dat["openmx"][38015];

# %%
locdat[1:20,1:4]

# %%
dt = alldat.dectree

# %%
dt[1][[9,5]]

# %%
dts = Dict(i=>sort(alldat.dectree[i], rev=true) for i in 1:5) 

# %%
collect(keys(dts[4]))

# %%
sickdays_by_age = Dict{Int,Array{Int,1}}()  # empty
fromconds_by_age = Dict{Int,Array{Int,1}}()  # empty
for i in 1:5
    sickdays_by_age[i] = [k[1] for k in collect(keys(dts[i]))]
    fromconds_by_age[i] = [k[2] for k in collect(keys(dts[i]))]
end
display(sickdays_by_age)
display(fromconds_by_age)

# %%
infected_idx = findall(locdat[:,1] .== 2);
sorted_by_sickday = sortperm(locdat[infected_idx,4], rev=true);

# %%
@elapsed for p in infected_idx[sorted_by_sickday]  # p is the index to the person
    (pstat, page, psickday, pcond) = locdat[p,[cpop_status, cpop_agegrp, cpop_sickday, cpop_cond]]

    sickdayfound = findall(x->x==psickday, sickdays_by_age[page])
    if isempty(sickdayfound)  # person's sickday doesn't match any decision point sickday
        # test against sickdaylim, then increment
        if psickday == sickdaylim
            @error "person made it to end of sickdaylim and was not removed"
        else
            locdat[p,cpop_sickday] += 1
        end
    else
        condfound = findall(x->x==pcond, fromconds_by_age[page][sickdayfound])
        if isempty(condfound)
            if psickday == sickdaylim
                @error "person made it to end of sickdaylim and was not removed"
            else
                locdat[p,cpop_sickday] += 1
            end
        else
            # do the transition for sickday and from cond and the probabilities of all outcomes at this branch
            dtkey = [psickday, pcond]
            probs = dts[page][dtkey]["probs"]
            outcomes = dts[page][dtkey]["outcomes"]
            choice = rand(Categorical(probs), 1)
            tocond = outcomes[choice][]
            if tocond in [dead, recovered]  # change status, leave cond and sickday as last state before death or recovery
                locdat[p, cpop_status] = tocond
            else   # change disease condition
                locdat[p, cpop_cond] = tocond
                locdat[p, cpop_sickday] += 1  
            end                      
        end
    end    
end

# %%
locdat[95515,1:4]
