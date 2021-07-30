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
using DataFrames
using Plots
using Printf
pyplot()

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :default, agegrps)

# %%
str_50 = sd_gen(start=50, comply=.8, cf=(.2,1.2), tf=(.18,.41))

# %%
# working with specific locale
seattle = 53033

# %%
alldat, series = run_a_sim(180, seattle, showr0=false, silent=true,
        runcases=[seed_1_6]);

# %%
cumplot(series,seattle,geo=alldat.geo)

# %%
sea_outcome = virus_outcome(series,seattle, base=:pop)
println(sea_outcome)
for k in keys(sea_outcome)
    @printf("%-12s %f\n", k, sea_outcome[k])
end

# %% [markdown]
# #### Death Percentage Across Age Groups

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
deadvals = series[seattle][:cum][end,[map2series.dead]...]
pctvals = round.([deadvals[i] / deadvals[6] for i in 1:length(deadvals)], digits=3)
death_dist_by_age = hcat(agelabels, deadvals, pctvals)

# %% [markdown]
# #### Death Percentage of Infected *Within* Each Age Group

# %%
dead = series[seattle][:cum][end, map2series.dead] 
infected = series[seattle][:cum][1,map2series.unexposed] .- series[seattle][:cum][end,map2series.unexposed]
death_pct_infected_within_age = round.(dead ./ infected, digits=5)
cats = hcat(agelabels, death_pct_infected_within_age)

# %% [markdown]
# #### Death Percentage of Population *Within* Each Age Group

# %%
pop = series[seattle][:cum][1,map2series.unexposed]
death_pct_bypop_within_age = round.(dead ./ pop, digits=5)
hcat(agelabels, death_pct_bypop_within_age)

# %% [markdown]
# #### Severe Percentage of Infected *Within* Each Age Group

# %%
sev_outcome = sum(clamp.(series[seattle][:new][:, map2series.severe], 0, 10_000_000), dims=1)'
sev_pct_infected_byage = round.(sev_outcome ./ infected, digits=5)
hcat(agelabels, sev_pct_infected_byage)

# %% [markdown]
# #### Severe Percentage of Population *Within* Each Age Group

# %%
sev_pct_pop_byage = round.(sev_outcome ./ pop, digits=5)
hcat(agelabels, sev_pct_pop_byage)

# %% [markdown]
# ### Recovered Distribution by Age Group

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
recovals = series[seattle][:cum][end,[map2series.recovered]...]
pctvals = round.([recovals[i] / recovals[6] for i in 1:length(recovals)], digits=3)
deadtbl = hcat(agelabels, recovals, pctvals)

# %% [markdown]
# ### Unexposed Percentage by Age Group

# %%
agelabels = ["0-20", "20-40", "40-60", "60-80", "80+", "Total"]
unexvals = series[seattle][:cum][end,[map2series.unexposed]...]
pctvals = round.([unexvals[i] / unexvals[6] for i in 1:length(unexvals)], digits=3)
deadtbl = hcat(agelabels, unexvals, pctvals)

# %% [markdown]
# ### Worldometers Death Demographics for New York City
#
# <img src=attachment:image.png width="500" height="500">

# %% [markdown]
# #### CDC Age Demographics for Covid-19 Deaths
# ##### through May 20, 2020 (based on slow reporting verified death reports--not latest)
# <img src="attachment:image.png" width="700px">

# %%
deadseries = series[seattle][:cum][:,[map2series.dead]...]
n = size(deadseries,1)

# %%
ageserieslabels = [agelabels[1] agelabels[2] agelabels[3] agelabels[4] agelabels[5]]
areaplot(1:n, deadseries[:,1:5],labels=ageserieslabels, title="Deaths by Age Group")

# %%
[deadseries[n,1:6] deadseries[n,1:6] ./ deadseries[n,6]]

# %% [markdown]
# ## Plots by Disease Condition

# %%
condseries = series[seattle][:cum][:,[map2series.nil[6], map2series.mild[6], map2series.sick[6], 
            map2series.severe[6]]]
n = size(condseries,1);

# %%
condlabels = ["nil", "mild", "sick", "severe"]
day = 180
condday = series[seattle][:cum][day,[map2series.nil[6], map2series.mild[6], map2series.sick[6], 
            map2series.severe[6]]]
condend = series[seattle][:cum][end,[map2series.nil[6], map2series.mild[6], map2series.sick[6], 
            map2series.severe[6]]]
condpct = round.(condday ./ sum(condday), digits=2)
println("Approximate Percentage Disease Condition\n(across all ages)")
condtbl = hcat(condlabels, condday, condpct)

# %%
condserieslabels = [condlabels[1] condlabels[2] condlabels[3] condlabels[4]]
areaplot(1:n, condseries[:,:],labels=condserieslabels, 
    title="Disease Conditions Over Time\n(across all ages)")

# %%
condserieslabels = [condlabels[4]]
areaplot(1:n, condseries[:,4],labels="Severe", title="Potential Hospital Burden")
maxsevere = maximum(condseries[:, 4])
half_yscale = floor(Int, maxsevere * 0.7)
annotate!((6,half_yscale,Plots.text("Burden: $maxsevere", 10, :left)))

# %%

# %%

# %%

# %% [markdown]
# ### Check the Basic Identities

# %%
cumhistmx = alldat.dat["cumhistmx"]
newhistmx = alldat.dat["newhistmx"]
openmx = alldat.dat["openmx"];

# %%
locale = seattle
outcome = (
           total_infected = series[locale][:cum][1, 6] - series[locale][:cum][180,6],
           total_pop = series[locale][:cum][180,6] + series[locale][:cum][180,54],
           whos_left = series[locale][:cum][180,map2series.dead[6]] + series[locale][:cum][180,map2series.recovered[6]]
              + series[locale][:cum][180,map2series.infectious[6]] + series[locale][:cum][180,map2series.unexposed[6]],
           end_unexposed = series[locale][:cum][180,map2series.unexposed[6]],
           end_infected = series[locale][:cum][180,map2series.infectious[6]],
           end_recovered = series[locale][:cum][180,map2series.recovered[6]],
           end_dead = series[locale][:cum][180,map2series.dead[6]]
       )

# %%
transeries = DataFrame(transq)


# %%
trans = (dead = sum(transeries[:,:dead]), recovered = sum(transeries[:,:recovered]))

# %%
err = outcome.total_infected - (trans.recovered + trans.dead + outcome.end_infected)

# %%
spreadseries = day2df(spreadq)
check_infected = sum(spreadseries[:,:infected])

# %% [markdown]
# end_exposed is ok (off by 2 from age rounding and 6 seeds)
# total infected is ok; matches check_infected

# %%

# %%

# %%

# %% [markdown]
# ##### Some random documentation for me

# %%
plotattr()

# %%
plotattr(:Subplot)

# %%
plotattr(:Series)

# %%
plotattr(:Plot)

# %%
plotattr("size")

# %%
plotattr(:Axis)

# %%
