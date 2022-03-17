# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.2
#   kernelspec:
#     display_name: Julia 1.7.0
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
export @ballocated, @belapsed, @benchmark, @benchmarkable, @benchmarkset, @bprofile, @btime, @case, @tagged, BenchmarkGroup, BenchmarkTools, addgroup!, allocs, gctime, improvements, invariants, isimprovement, isinvariant, isregression, judge, leaves, loadparams!, mean, median, memory, params, ratio, regressions, rmskew, rmskew!, trim, tune!, warmup
using Distributions
using YAML
using PrettyPrint
using OrderedCollections
using Revise
Revise.revise()
using LinearAlgebra

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# # Test setup and population matrix

# %% tags=[]
ndays = 180
locale = 38015
model = buildsim(ndays, locale;  
    dovax = true,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantfilename = "variants.yml",
);

# %%
keys(model)

# %% [markdown]
# ### Data tables

# %%
locdat = model.dat["popdat"][locale]

# %%
ages = model.dat["agegrp_idx"][locale]

# %%
columnnames(locdat)

# %%
countmap(locdat.agegrp)

# %%
countmap(locdat.status)  # everyone begins as unexposed

# %% [markdown]
# ### Social Parameters

# %%
fieldnames(typeof(model.social))

# %%
typeof(model.social.gammashape)

# %%
model.social.contactfactors

# %%
model.social.contactfactors[Int(sick)-4, Int(age0_19)]

# %%
touchfactors = model.social.touchfactors

# %% [markdown]
# ### Parameters for infection based on variant

# %%
fieldnames(typeof(model.infectset[:omicron]))

# %%
model.infectset

# %%
pprintln(model.infectset)

# %%
pprintln(model.infectset[:base])

# %%
sendbase = model.infectset[:base].sendrisk
recvbase = model.infectset[:base].recvrisk

# %%
newomisend = zeros(25)
newomisend[1:11] = 1.5 .* sendbase[1:11]
newomisend

# %%
pprint(model.infectset[:alpha])

# %%
pprintln(model.infectset[:omicron])

# %%
# is shifter working?
shifter(touchfactors, (.18, .3)...)[:, Int(age40_59)]

# %%
model.transitionset # the decision transition matrices for all age groups are loaded

# %%
basetransition = model.transitionset[:base]

# %%
fieldnames(typeof(basetransition))

# %%
fieldnames(typeof(basetransition.tree))

# %%
age80tree = basetransition.tree.age80_up[1].sickday

# %%
age80tree = basetransition.tree.age80_up[1].transition

# %%
keys(age80tree)  # array of Transitiondef

# %%
basetransition.factors

# %%
basetransition.factors.vaxhalflifeadjust

# %%
println(typeof(age80tree))
breakday_idx = 5
println("for breakday $breakday_idx")
println("fieldnames: ", fieldnames(typeof(age80tree[breakday_idx])))
println("sickday ",age80tree[breakday_idx].sickday)
println("transition \n", age80tree[breakday_idx].transition)

# %% [markdown]
# ### Sanity check the transition tree

# %%
cs.sanitycheck(basetransition.tree)

# %% [markdown]
# # Run a simulation

# %% [markdown]
# ### Build the simulation model

# %%
ndays = 720
locale = 38015
model = buildsim(ndays, locale;  
    dovax = true,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantfilename = "variants.yml",
);

# %% [markdown]
# ### Create a seed case

# %%
seed20_39_day1 = makesickseedfunc(; cond=nil, variant=:base, sickday=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=1, forstartofday=true)
seed40_59_day1 = makesickseedfunc(; cond=nil, variant=:base, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=1, forstartofday=true)                            

# %%
seed20_39_omicron = makesickseedfunc(; cond=nil, variant=:omicron, sickday=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=360, forstartofday=true)
seed40_59_omicron = makesickseedfunc(; cond=nil, variant=:omicron, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=360, forstartofday=true)                            

# %% [markdown]
# ### Run the simulation model

# %%
popdat, series = runsim(model;
            dovax=true,
            dovariant = true,
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_omicron, seed40_59_omicron]   # or seed_1_6 if using old way
            );

# %% [markdown]
# ### Plot results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orange line labeled Infectious, which shows the current number of infected people, is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: There were some sick people as of the day before. Some more people got sick today. Some people got better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. So net infected is yesterday + new today - recovered today - died today. Newspaper tracking shows the new infections of each day--who got sick today? Tomorrow, if no one new got sick the line would be at zero--even though the people who got sick yesterday aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. 

# %%
cumplot(series, locale, [:nil, :mild, :sick, :severe])

# %%
newplot(series, locale)

# %%
cumplot(series, locale, [:base, :omicron])

# %%
colmap = series[locale].cols

# %%
series[locale].cum[end, colmap[:base]]

# %%
series[locale].cum[end, colmap[:base]]

# %%
(cs.vaxlist, cs.variantlist)

# %%
keys(model)

# %%
keys(popdat)

# %%
keys(series)

# %%
model.dat

# %%
popdat[locale]

# %%
locdat = popdat[locale]
countmap(locdat.cond)

# %%
countmap(locdat.status)

# %%
countmap(locdat.variant)

# %%
virus_outcome(series, locale, base=:pop)

# %%
fieldnames(typeof(series[locale]))

# %%
map2series = series[locale].cols

# %%
println(":cum ", typeof(series[locale].cum))
println(":cum ", size(series[locale].cum))
println(":new ", size(series[locale].new))

# %%
series[locale].new[:, map2series[:infectious]]

# %% [markdown]
# ## Test a social distancing case

# %%
sd1 = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    

# %%
sd1_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, include_ages=[])

# %%
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, runcases=[seed_1_6, sd1, sd1_end]);

# %%
virus_outcome(series, locale, base=:pop)

# %%
cumplot(series, locale)

# %%
cumplot(series, locale,[:infectious, :dead])

# %%
outdat = result_dict.dat["popdat"][locale]
all(outdat.sdcomply .== :none)

# %% [markdown]
# ## Social distancing only among those age40_59, age60_79, age80_plus

# %%
sdolder = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# %%
sdolder_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# %%
popdat, series = runsim(model; 
                    runcases=[seed_1_6, sdolder, sdolder_end]);


# %%
olderdat =popdat[locale]

sd = findall(olderdat.sdcomply .!= :none)

@Select(agegrp, sdcomply)(olderdat[sd])

# %%
cumplot(series, locale)

# %%
cumplot(series, locale, [:infectious, :dead])

# %% [markdown]
# ## Social Distancing starts with everyone and then the younger folks party

# %%
sdyoung_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age0_19, age20_39])    

# %%
popdat, series = runsim(model,
    runcases=[seed_1_6, sd1, sdyoung_end]);

# %%
mixdat = popdat[locale]

sd = findall(mixdat.sdcomply .!= :none)
young = findall((mixdat.agegrp .== age0_19) .| (mixdat.agegrp .== age20_39))
sd_young_idx = intersect(sd, young)
old = findall((mixdat.agegrp .== age40_59) .| (mixdat.agegrp .== age60_79) .| (mixdat.agegrp .== age80_up));

# %%
youngtab = Table(mixdat[young])
count(youngtab.sdcomply .== :none)

# %%
oldtab = Table(mixdat[old])
count(oldtab.sdcomply .!= :none)

# %%
typeof(mixdat.agegrp)

# %%
mixdat = popdat[locale]

# %%
cumplot(series, locale, [:infectious, :dead])

# %%
@Select(status, agegrp, cond, sdcomply)(locdat)
