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
using CovidSim

# %%
seattle = (;fips=53033); newyork=(;fips=36061); bismarck=(;fips=38015)

# %%
geo = CovidSim.readgeodata("../data/geo2data.csv")
geo[:,1:7]

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :base, agegrps)

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       runcases=[seed_1_6]);

# %%
nycfull=cumplot(series,newyork.fips,geo=geo)

# %%
infection_outcome(series, newyork.fips)

# %% [markdown]
# ### Strong Social Distancing

# %% [markdown]
# Reset the model to defaults.

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       runcases=[seed_1_6]);

# %%
str_50 = sd_gen(start=50, comply=.9, cf=(.5,1.2), tf=(.18,.42))

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50],
       runcases=[seed_1_6]);

# %%
nycsd=cumplot(series,newyork.fips,[recovered, infectious, dead],geo=geo)

# %%
infection_outcome(series, newyork.fips)

# %% [markdown]
# ### Open Totally (which won't happen...)
# This uses opening back to essentially no social distancing and an R0 between 1.6 and 1.8. People will voluntarily be more prudent and government recommendations and policies will provide for more limited opening. So, this shows why complete opening isn't possible:  the full force of the virus does return with only a slight delay.

# %%
# Reset to defaults
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
open = sd_gen(start=95, comply=0.0, cf=(.3,1.8), tf=(.18,.62))

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50, open],
       runcases=[seed_1_6]);

# %%
nycopen=cumplot(series,newyork.fips,[infectious, dead],geo=geo)

# %%
r0_sim(;sa_pct=[1.0,0.0,0.0], density_factor=1.25, dt=alldat.dectree, cf=[], tf=[],
                compliance=[1.0], shift_contact=(.6,1.8), shift_touch=(.18,.62), disp=false, spreadparams=spreadparams)

# %% [markdown]
# The preceding estimate of R0 tracks one cohort across all demographic groups of the simulation through the duration of the disease for 25 days, though much of the cohort will "drop out" through recovery or death prior to the 25th day. Infectiveness varies across the number of days a person has been infected.  

# %% [markdown]
# ### Limited Opening

# %%
# reset the model to defaults
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[],
       runcases=[seed_1_6]);

# %%
open = sd_gen(start=95, comply=0.7, cf=(.5,1.5), tf=(.25,.50))

# %%
alldat, series = run_a_sim(180, newyork.fips, showr0=true, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50, open],
       runcases=[seed_1_6]);

# %%
nycopen2=cumplot(series,newyork.fips,geo=geo)

# %% [markdown]
# Even with some restrictions still in place, this opening results in only about 10 to 11% reduction in deaths compared to never having implemented social distancing or social distancing followed by completely opening. This alternative would *not* be recommended. Rigorous implementation of test, trace and isolate with idespread testing of asymptomatic people, same-day test results, and high compliance with contact tracing and quarantine isolation is very hard to achieve, but is a preferred alternative.  If the high degree of compliance and implementation rigor is not possible, is there an alternative?

# %% [markdown]
# # An Alternative: Fewer Restrictions with Protection of the Vulnerable
# Note that this cases models vulnerable people as the age groups from 60 to 80 and over 80. Other people are vulnerable: people with diabetes, hypertension, immuno-compromise, cancer patients, smokers, and others across age groups. It's beyond this model to attempt to represent these vulnerabilities across age groups at this poing.

# %%
open_more = sd_gen(start=95, cf=(.5,1.55), tf=(.25,.55),comply=.6)

# %%
r0_sim(;sa_pct=[1.0,0.0,0.0], density_factor=1.25, dt=alldat.dectree, cf=[], tf=[],
                compliance=[.65], shift_contact=(.5,1.55), shift_touch=(.25,.52), disp=false, spreadparams=spreadparams)

# %%
function isolate_vulnerable(locale; opendat=openmx, isodat=isolatedmx,testdat=openmx, spreadparams=spreadparams)
    if day_ctr[:day] == 105
        isolate!(.70,[unexposed, nil,mild,sick, severe],[5],1:durationlim, locale; opendat=opendat, isodat=isodat)
        isolate!(.50,[unexposed,nil,mild,sick, severe],[4],1:durationlim, locale; opendat=opendat, isodat=isodat)
    end
end

# %%
# reset the model to defaults
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       runcases=[]);

# %%
alldat, series = run_a_sim(180, newyork.fips, showr0=true, silent=true,
       dtfilename="../parameters/dec_tree_all_25.csv",
       spreadcases=[str_50, open],
       runcases=[seed_1_6, isolate_vulnerable]);

# %%
cumplot(series,newyork.fips,geo=geo)

# %% [markdown]
# **Assessment:**
# With opening somewhat more liberal and with lower compliance, we observe simulated deaths that are less than half of both the case of no social distancing ever and social distancing followed by opening with light restrictions held in place. The number of deaths is roughly 30 to 35% higher than maintaining restrictive social distancing for four months, until the end of the simulation duration of 180 days. But, that is not occurring anywhere in the US and is not politically possible nor practicable as compliance would continue to wain despite the terrible outcomes of opening up. 
#
# Comparing the protect the vulnerable approach to versions of opening up, with and without test, trace and isolate
# - deaths roughly 50% lightly restricted opening up;
# - somewhat lighter restrictions than the preceding;
# - deaths roughly 10% higher than very rigorous test, trace and isolate
# - deaths roughly 50% lower than loosely implemented test, trace and isolate
#

# %%
protected_group = age_dist .* series[newyork.fips][:cum][1,map2series.unexposed[6]]
protected_group = protected_group[4] + protected_group[5]

# %%
caseload_protected = 150 # per month
workers_protected = protected_group / caseload_protected
tests_required_month = 8 * workers_protected
tests_required_day = tests_required_month / 30

# %%
tests_required_month = 30 * workers_protected
tests_required_day = tests_required_month / 30

# %%
