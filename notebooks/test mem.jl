# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML
using PrettyPrint
using LinearAlgebra
using Dates
using Plots
using TableView
#using PlotlyBase


# %%
using CovidSim_ilm


# %%
cs = CovidSim_ilm


# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/src"))

# %% [markdown]
# # Test setup and population matrix
# %% [markdown]
# # Run a simulation
# %% [markdown]
# ### Build the simulation model

# %%
ndays = 850
locale = 38015
model = buildsim(ndays, locale;  
    day1 = Date("2020-01-01", "yyyy-mm-dd"),
    dovax = true,
    paramdir = "../sample_parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantfilename = "variants.yml",
);


# %%
seed20_39_day1 = makesickseedfunc(; cond=nil, variant=:base, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true);
seed40_59_day1 = makesickseedfunc(; cond=nil, variant=:base, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true);           


# %%
seed20_39_delta = makesickseedfunc(; cond=nil, variant=:delta, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
        cnt=3, forlocale=0, triggerdate=300, forstartofday=true);
seed40_59_delta = makesickseedfunc(; cond=nil, variant=:delta, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
        cnt=3, forlocale=0, triggerdate=300, forstartofday=true);        


# %%
seed20_39_omicron = makesickseedfunc(; cond=nil, variant=:omicron_ba1, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=660, forstartofday=true);
seed40_59_omicron = makesickseedfunc(; cond=nil, variant=:omicron_ba1, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=660, forstartofday=true);                            


# %%
seed20_39_omicron_ba2 = makesickseedfunc(; cond=nil, variant=:omicron_ba2, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=6, forlocale=0, triggerdate=690, forstartofday=true);
seed40_59_omicron_ba2 = makesickseedfunc(; cond=nil, variant=:omicron_ba2, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=6, forlocale=0, triggerdate=690, forstartofday=true);                            

# %% [markdown]
# ### Run the simulation model

# %%
popdat, series = runsim(model;
            dovax=true, vaxscheds=:loc38015,
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, seed20_39_omicron, seed40_59_omicron,
                      seed20_39_omicron_ba2, seed40_59_omicron_ba2]   # or seed_1_6 if using old way
            );
locdat = popdat[locale];


# %%
TableView.showtable(locdat[95000:95200])


# %%
TableView.showtable(series[locale].cum)


# %%
stat1 = cs.stat1(series, locale)


# %%
size_popdat = Base.summarysize(popdat)
size_series = Base.summarysize(series)
@show size_popdat size_series

# %% [markdown]
# ### Plot results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orange line labeled Infectious, which shows the current number of infected people, is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: There were some sick people as of the day before. Some more people got sick today. Some people got better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. So net infected is yesterday + new today - recovered today - died today. Newspaper tracking shows the new infections of each day--who got sick today? Tomorrow, if no one new got sick the line would be at zero--even though the people who got sick yesterday aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. 

# %%
cumplot(series, locale, [:nil, :mild, :sick, :severe])


# %%
newplot(series, locale, [:dead])


# %%
# png(joinpath(homedir(),"Desktop", "daily change.png"))


# %%
cumplot(series, locale, [:base, :delta, :omicron_ba1, :omicron_ba2])

# %% [markdown]
# ## Examine Vaccination process and outcomes

# %%
vax = :Pfizer
pfdose2 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vax]) .== 2)
pfdose1 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vax]) .== 1)
println("Count people with $vax 1 dose: $pfdose1 2 doses: $pfdose2 Any: $(pfdose1 + pfdose2)")
println("Doses used = $(2 * pfdose2 + pfdose1)")


# %%
vax = :Moderna
modose2 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vax]) .== 2)
modose1 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vax]) .== 1)
println("Count people with $vax 1 dose: $modose1 2 doses: $modose2 Any: $(modose1 + modose2)")
println("Doses used = $(2 * modose2 + modose1)")


# %%
vax = :JnJ
# jjdose2 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vax]) .== 2)
jjdose1 = count(length.(locdat.vaxrcvd[last.(locdat.vaxrcvd) .== vax]) .== 1)
println("Count people with $vax 1 dose: $jjdose1 Any: $(jjdose1)")
println("Doses used = $(jjdose1)")


# %%
println("People with any vaccine: $(pfdose2+pfdose1+modose2+modose1+jjdose1)")


# %%
newplot(series, locale, [:Pfizer], days=380:700)


# %%
newplot(series, locale, [:Moderna], days=380:700)


# %%
newplot(series, locale, :JnJ, days=380:700)


# %%
newplot(series, locale, :unexposed)



# %% [markdown]
# ## How many people got sick multiple times?
# 
# 

# %%
sort(countmap(length.(locdat.variant)))


# %%
@Select(vaxday, sickday)(locdat[(last.(locdat.sickday) .> 0) .& (first.(locdat.vaxday).> 0) .& (last.(locdat.sickday) .> first.(locdat.vaxday))], )

# %% [markdown]
# ## Test a social distancing case

# %%
sd1 = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    


# %%
sd1_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, include_ages=[])


# %%
popdat, series = runsim(model;
            dovax=true, vaxscheds=:loc38015, showr0=false, silent=true, 
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, seed20_39_omicron, seed40_59_omicron,
                      seed20_39_omicron_ba2, seed40_59_omicron_ba2, sd1, sd1_end]);
locdat = popdat[locale];


# %%
vo=virus_outcome(series, locale, base=:pop)
display(vo)


# %%
sum(values(vo))


# %%
cumplot(series, locale)


# %%
cumplot(series, locale,[:infectious, :dead])


# %%
cumplot(series, locale,[:base, :delta, :omicron_ba1, :omicron_ba2])

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
                        dovax=true, vaxscheds=:loc38015,
                        runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, 
                                  seed20_39_omicron, seed40_59_omicron,
                                  seed20_39_omicron_ba2, seed40_59_omicron_ba2, 
                                  sdolder, sdolder_end]);
locdat = popdat[locale];


# %%
olderdat =popdat[locale]

sd = findall(olderdat.sdcomply .!= :none)

@Select(agegrp, sdcomply)(olderdat[sd])
# @Select(agegrp, sdcomply)(olderdat)


# %%
cumplot(series, locale)


# %%
freemem_plot3 = (Sys.free_memory() / 2^20)


# %%
cumplot(series, locale, [:infectious, :dead])

# %% [markdown]
# ## Social Distancing starts with everyone and then the younger folks party

# %%
sdyoung_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age0_19, age20_39])    


# %%
popdat, series = runsim(model,
                        dovax=true, vaxscheds=:loc38015,
                        runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, 
                                  seed20_39_omicron, seed40_59_omicron,
                                  seed20_39_omicron_ba2, seed40_59_omicron_ba2, 
                                  sd1, sdyoung_end]);
locdat = popdat[locale];


# %%
freemem_run4 = (Sys.free_memory() / 2^20)


# %%
@show freemem_start freemem_pkg freemem_build freemem_run1 freemem_plot1 freemem_run2 freemem_plot2 freemem_run3 freemem_plot3 freemem_run4


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


