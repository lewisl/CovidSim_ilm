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
cd(joinpath(homedir(),"Dropbox/Online Coursework/Covid/ilm-src")); 

# %%
using Plots
pyplot()

# %%
using CovidSim_ilm

# %%
bismarck = (;fips=38015)
newyork = (;fips=36061)

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :default, agegrps)

# %%
alldat, series = run_a_sim(180, bismarck.fips, showr0=false, silent=true,
        spreadcases=[],
        runcases=[seed_1_6]);
geo = alldat.geo;

# %% [markdown]
# # Smaller City: No social distancing

# %%
bisfull=cumplot(series,bismarck.fips,geo=geo)

# %%
virus_outcome(series, bismarck.fips)

# %%
str_50 = sd_gen(start=50, comply=.9, cf=(.5,1.2), tf=(.18,.42))
alldat, series = run_a_sim(180, bismarck.fips, showr0=false, silent=true,
    spreadcases=[str_50],
    runcases=[seed_1_6]);

# %% [markdown]
# # Moderately Strong Social Distancing starts on Day 50

# %%
bissd50=cumplot(series,bismarck.fips,geo=geo)

# %%
virus_outcome(series, bismarck.fips)

# %%
open = sd_gen(start=80, comply=0.7, cf=(.5,1.5), tf=(.25,.50))

# %%
alldat, series = run_a_sim(180,bismarck.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6]);

# %% [markdown]
# # With low cases and deaths very low-->open up significantly, but not completely

# %%
bisopen=cumplot(series,bismarck.fips,geo=geo)

# %%
bisbump=cumplot(series,bismarck.fips,[infectious, dead],geo=geo)

# %%
virus_outcome(series, bismarck.fips)

# %% [markdown]
# # Adopt Test, Trace and Isolate

# %%
t_n_t100_160=CovidSim_ilm.t_n_t_case_gen(100,160,tc_perday=1000, test_delay=3,
    breakout_pct=0.2, q_comply=0.75,past_contacts=false)

# %%
alldat, series = run_a_sim(180, bismarck.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %% [markdown]
# # Apply Test, Trace and Isolate 
# - test capacity of 1% of population per day
# - quarantine compliance of 75%
# - delay of 3 days returning test results

# %%
cumplot(series,38015,geo=geo)

# %% [markdown]
# # Apply Test, Trace and Isolate--raising test capacity to 5% of population
#
# - 4800 tests per day (roughly 5% of population)
# - quarantine compliance of 95%
# - test results returned same day

# %%
t_n_t100_160=CovidSim_ilm.t_n_t_case_gen(100,160,tc_perday=4800, test_delay=0,
    breakout_pct=0.20, q_comply=0.95, past_contacts=false)

# %%
alldat, series = run_a_sim(180,bismarck.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %%
cumplot(series,bismarck.fips,geo=geo)

# %% [markdown]
# # Apply Test, Trace and Isolate
# - testing capacity of 10% of the population
# - quarantine compliance of 97%
# - delay in obtaining test results of 0 days

# %%
t_n_t100_160=CovidSim_ilm.t_n_t_case_gen(100,160,tc_perday=10_000, test_delay=0,
    breakout_pct=0.20, q_comply=0.97,past_contacts=false)

# %%
alldat, series = run_a_sim(180, bismarck.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %%
cumplot(series,bismarck.fips,geo=geo)

# %%

# %%

# %% [markdown]
# # Very large city: no social distancing

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
    spreadcases=[],
    runcases=[seed_1_6]);

# %%
nycfull=cumplot(series, newyork.fips, geo=geo)

# %%
infection_outcome(series, newyork.fips)

# %%
str_50 = sd_gen(start=50, comply=.9, cf=(.5,1.2), tf=(.18,.42))
alldat, series = run_a_sim(180, newyork.fips, showr0=false, silent=true,
    spreadcases=[str_50],
    runcases=[seed_1_6]);

# %% [markdown]
# # Moderately Strong Social Distancing starts on Day 50

# %%
nycsd50=cumplot(series,newyork.fips,geo=geo)

# %% [markdown]
# Note: The total of deaths is a bit higher than NYC deaths as of May 23 at 21086, which is day 122, while the simulation is run for almost another 2 months.

# %%
virus_outcome(series, newyork.fips)

# %%
open = sd_gen(start=95, comply=0.7, cf=(.5,1.5), tf=(.25,.50))

# %%
alldat, series = run_a_sim(180, newyork.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6]);

# %% [markdown]
# # With declining cases and declining deaths-->open up later on day 95

# %%
nycopen=cumplot(series,newyork.fips,geo=geo)

# %%
nycbump=cumplot(series, newyork.fips, [infectious, dead],geo=geo)

# %% [markdown]
# Note: We don't see a "double-dip" pattern for New York City because the original curve of infection had not turned as sharply down, but had moderated considerably. When opening occurs, the infection rate looks like unrestrained virus transmission.

# %%

# %%
t_n_t100_160=CovidSim_ilm.t_n_t_case_gen(100,160,tc_perday=80_000,breakout_pct=0.20, q_comply=0.75,past_contacts=false)

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %% [markdown]
# # Apply Test, Trace and Isolate: 
# - 80,000 tests per day 
# - quarantine compliance of 75%
# - delay in returning test results of 3 days

# %%
nyc_tnt_min = cumplot(series,newyork.fips,geo=geo)

# %% [markdown]
# # Apply Test, Trace and Isolate 
#
# - 835,000 tests per day (roughly 10% of population)
# - quarantine compliance of 95%
# - test results returned same day

# %%
t_n_t100_160=CovidSim_ilm.t_n_t_case_gen(100,160,tc_perday=835_000, test_delay=0,
    breakout_pct=0.20, q_comply=0.95, past_contacts=false)

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %%
nyc_tnt_max = cumplot(series,newyork.fips,geo=geo)

# %% [markdown]
# Note: returning tests results same day results in 25% reduction in deaths compared to a delay of 3 days. 

# %% [markdown]
# # Apply Test, Trace and Isolate--raising test capacity to 5% of population
#
# - 420,000 tests per day (roughly 5% of population)
# - quarantine compliance of 95%
# - test results returned same day

# %% jupyter={"source_hidden": true}
t_n_t100_160=CovidSim_ilm.t_n_t_case_gen(100,160,tc_perday=420_000, test_delay=0,
    breakout_pct=0.20, q_comply=0.95, past_contacts=false)

# %%
alldat, series = run_a_sim(180,newyork.fips, showr0=false, silent=true,
    spreadcases=[str_50, open],
    runcases=[seed_1_6, t_n_t100_160]);

# %%
nyc_tnt_mid = cumplot(series,newyork.fips,geo=geo)

# %%
