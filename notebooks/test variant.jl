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
model180 = buildsim(ndays, locale;  
    dovax = true,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantfilename = "variants.yml",
);

# %%
keys(model180)

# %% [markdown]
# ### Data tables

# %%
locdat = model180.dat["popdat"][locale]

# %%
ages = model180.dat["agegrp_idx"][locale]

# %%
columnnames(locdat)

# %%
countmap(locdat.agegrp)

# %%
countmap(locdat.status)  # everyone begins as unexposed

# %% [markdown]
# ### Series

# %%
fieldnames(typeof(model180.series))

# %%
model180.series[locale].cum

# %%
for col in columnnames(model180.series[locale].cum)
    println(col)
end

# %% [markdown]
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
model180.infectset

# %%
fieldnames(typeof(model180.infectset[:omicron_ba1]))

# %%
CovidSim_ilm.variantlist

# %%
pprintln(model180.infectset)

# %%
pprintln(model180.infectset[:base])

# %%
sendbase = model180.infectset[:base].sendrisk
recvbase = model180.infectset[:base].recvrisk

# %%
newomisend = zeros(25)
newomisend[1:11] = 1.5 .* sendbase[1:11]
newomisend

# %%
pprint(model180.infectset[:alpha])

# %%
pprintln(model180.infectset[:omicron_ba1])

# %%
# is shifter working?
shifter(touchfactors, (.18, .3)...)[:, Int(age40_59)]

# %%
model180.transitionset # the decision transition matrices for all age groups are loaded

# %%
basetransition = model180.transitionset[:base]

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
# ## Vaccines and Vaccination Schedule

# %%
pprintln(model180.vaxset)

# %%
pprintln(model180.vaxschedset)

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
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true)
seed40_59_day1 = makesickseedfunc(; cond=nil, variant=:base, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true)                            

# %%
seed20_39_omicron = makesickseedfunc(; cond=nil, variant=:omicron_ba1, sickday=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=360, forstartofday=true)
seed40_59_omicron = makesickseedfunc(; cond=nil, variant=:omicron_ba1, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=360, forstartofday=true)                            

# %%
seed20_39_omicron_ba2 = makesickseedfunc(; cond=nil, variant=:omicron_ba2, sickday=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=460, forstartofday=true)
seed40_59_omicron_ba2 = makesickseedfunc(; cond=nil, variant=:omicron_ba2, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=460, forstartofday=true)                            

# %% [markdown]
# ### Run the simulation model

# %%
popdat, series = runsim(model;
            dovax=true,
            dovariant = true,
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_omicron, seed40_59_omicron,
                      seed20_39_omicron_ba2, seed40_59_omicron_ba2]   # or seed_1_6 if using old way
            );
locdat = popdat[locale];

# %%
filt_vaccinated = findall(last.(popdat[38015].vaxrcvd) .!= :none)
if length(filt_vaccinated) > 0
    vax_today = @inbounds countmap(last.(popdat[38015].vaxrcvd[filt_vaccinated]))  # keys are symbol
else
    vax_today = Dict()
end
println(vax_today)
println("Doses given: ",sum(values(vax_today)))


# %%
vaxsym = :Moderna
dose1 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vaxsym]) .== 1)
dose2 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vaxsym]) .== 2)
println(vaxsym, " 1 dose cnt: ", dose1, " 2 dose cnt: ", dose2, " total ", 2*dose2 + dose1, " people cnt: ", dose1 + dose2)

# %% [markdown]
# ### Plot results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orange line labeled Infectious, which shows the current number of infected people, is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: There were some sick people as of the day before. Some more people got sick today. Some people got better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. So net infected is yesterday + new today - recovered today - died today. Newspaper tracking shows the new infections of each day--who got sick today? Tomorrow, if no one new got sick the line would be at zero--even though the people who got sick yesterday aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. 

# %%
cumplot(series, locale, [:Pfizer, :Moderna, :JnJ, :totvaccinated], days=160:300)

# %%
cumplot(series, locale, [:nil, :mild, :sick, :severe])

# %%
newplot(series, locale, [:unexposed])

# %%
cumplot(series, locale, [:base, :omicron_ba1, :omicron_ba2])

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
locdat = popdat[locale]
@btime countmap($locdat.cond)

# %%
@btime cnt_cond($locdat.cond)

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
