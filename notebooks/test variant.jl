# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Julia (4 threads) 1.8.1
#     language: julia
#     name: julia-4-threads-1.8
# ---

# %%
using Pkg
# Pkg.activate("/Users/lewis/.julia/environments/v1.8")

# %%
readdir()

# %% tags=[]
# Run this when actively modifying the code in the package
Pkg.develop(path="/Users/lewis/Dropbox/Covid Modeling/Covid-ILM")

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
using Plots
using TableView
using SplitApplyCombine
using Tables
using PrettyTables

# %%
# Use this command running locally
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/src"))

# %%
# Use this command running on Github Codespaces
cd("/workspaces/CovidSim_ilm/src")

# %%
;pwd

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[]
# # Test setup and population matrix

# %% tags=[]
ndays = 180
locale = 38015
model180 = buildsim(ndays, locale;  
    day1 = Date("2020-01-01", "yyyy-mm-dd"),
    dovax = true,
    paramdir = "../sample_parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    scheddir = "vaccine_100k",
    variantfilename = "variants.yml",
);

# %% [markdown]
# # Run the small model

# %% [markdown]
# ### Seed with base variant on day 1

# %%

seed20_39_day1 = makesickseedfunc(; cond=nil, variant=:base, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true);
seed40_59_day1 = makesickseedfunc(; cond=nil, variant=:base, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true);           


# %% [markdown] tags=[]
# ### Run and Compile

# %%
popdat, series = runsim(model180;
            dovax=false, vaxscheds=:all,
            runcases=[seed20_39_day1, seed40_59_day1]   
            );
locdat = popdat[locale];

# %% [markdown] tags=[]
# ### Plot results

# %%
cumplot(series, locale, [:unexposed, :infectious, :recovered, :dead])

# %% [markdown] tags=[]
# ### Run and profile
#
# # only in local VS Code
# @profview popdat, series = runsim(model180;
#                     dovax=false, vaxscheds=:all,
#                 runcases=[seed20_39_day1, seed40_59_day1]   
#                 );

# %% [markdown]
# ### R0 Simulation

# %%
r0 = cs.r0_sim(200_000, AGE_DIST, model180.progressionset, model180.trvec, model180.infectset, 
                model180.vaxset, model180.social, :base, 1.0, 3)

# %%
keys(model180)

# %%
typeof(model180.seriescolnames)

# %%
model180.seriescolnames

# %%
scn = model180.seriescolnames
@btime scn[:condcols][Symbol(sick)][Symbol(age0_19)]


# %%
@btime Symbol(:sick, "_", age0_19)

# %% [markdown]
# ### Data tables

# %%
locdat = model180.dat.popdat[locale]

# %%
ages = model180.dat.agegrp_idx[locale]

# %%
columnnames(locdat)

# %%
countmap(locdat.agegrp)

# %%
countmap(locdat.status)  # everyone begins as unexposed

# %% [markdown]
# ### Series

# %%
typeof(model180.series) <: Dict

# %%
lochist = model180.series[locale];
newhist = lochist.new;
cumhist = lochist.cum;

# %%
newhist

# %%
for col in columnnames(lochist.cum)
    println(col)
end

# %% tags=[]
age = "total"
symb_age = Symbol(age)
symb_status = Symbol(unexposed)
cumcoll = getproperty(lochist.cum, Symbol(symb_status, "_", symb_age))

# %%
thisday = 100
cumcell = getproperty(lochist.cum, Symbol(symb_status, "_", symb_age))[thisday]

# %% [markdown]
# ### Geo Data

# %%
model180.geo

# %%
typeof(model180.geo) <: Table

# %% [markdown] tags=[]
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

# %% [markdown]
# ### Parameters for infection based on variant

# %%
model180.infectset

# %%
fieldnames(typeof(model180.infectset[:omicron_ba1]))

# %%
pprintln(model180.infectset)

# %%
pprintln(model180.infectset[:base])

# %%
sendbase = model180.infectset[:base].sendrisk
recvbase = model180.infectset[:base].recvrisk

# %%
pprint(model180.infectset[:alpha])

# %%
pprintln(model180.infectset[:omicron_ba1])

# %%
# is shifter working?
shifter(touchfactors, (.18, .3)...)[:, Int(age40_59)]

# %% [markdown]
# ## Transitionset

# %%
model180.transitionset # the decision transition matrices for all age groups are loaded

# %%
basetransition = model180.transitionset[:base]
fieldnames(typeof(basetransition))

# %%
fieldnames(typeof(basetransition.tree))

# %%
display(basetransition.tree.age0_19)

# %%
cs.sanitycheck(basetransition.tree)

# %%
alphatransition = model180.transitionset[:alpha].tree

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

# %%

# %% [markdown]
# ## Vaccines and Vaccination Schedule

# %%
model180.vaxset

# %%
model180.vaxschedset[:loc38015]

# %% [markdown]
# ### Sanity check the transition tree

# %%
cs.sanitycheck(basetransition.tree)

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
    scheddir = "vaccine_100k",
    variantfilename = "variants.yml",
);

# %% [markdown]
# ### Seed with base variant on day 1

# %%

seed20_39_day1 = makesickseedfunc(; cond=nil, variant=:base, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true);
seed40_59_day1 = makesickseedfunc(; cond=nil, variant=:base, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=1, forstartofday=true);           

# %% [markdown]
# ### Seed with delta variants on day 300

# %%

seed20_39_delta = makesickseedfunc(; cond=nil, variant=:delta, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
        cnt=3, forlocale=0, triggerdate=300, forstartofday=true);
seed40_59_delta = makesickseedfunc(; cond=nil, variant=:delta, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
        cnt=3, forlocale=0, triggerdate=300, forstartofday=true);        

# %% [markdown]
# ### Seed with omicron_ba1 variant on day 660

# %%

seed20_39_omicron = makesickseedfunc(; cond=nil, variant=:omicron_ba1, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=640, forstartofday=true);
seed40_59_omicron = makesickseedfunc(; cond=nil, variant=:omicron_ba1, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, triggerdate=640, forstartofday=true);                            

# %% [markdown]
# ### Seed with omicron_ba2 on day 690

# %%

seed20_39_omicron_ba2 = makesickseedfunc(; cond=nil, variant=:omicron_ba2, duration=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=6, forlocale=0, triggerdate=680, forstartofday=true);
seed40_59_omicron_ba2 = makesickseedfunc(; cond=nil, variant=:omicron_ba2, duration=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=6, forlocale=0, triggerdate=680, forstartofday=true);                            

# %% [markdown]
# ### Seed with omicron_ba4_5 on day 740

# %%

seed20_39_omicron_ba4_5 = makesickseedfunc(; cond=nil, variant=:omicron_ba4_5, duration=1, filter=[ Term(:agegrp, age20_39), Term(:status, recovered) ], 
                            cnt=6, forlocale=0, triggerdate=740, forstartofday=true);
seed40_59_omicron_ba4_5 = makesickseedfunc(; cond=nil, variant=:omicron_ba4_5, duration=1, filter=[ Term(:agegrp, age40_59), Term(:status, recovered) ], 
                            cnt=6, forlocale=0, triggerdate=740, forstartofday=true);                            

# %% [markdown]
# ### Run the simulation model

# %%
popdat, series = runsim(model;
            dovax=true, vaxscheds=:all,
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, seed20_39_omicron, seed40_59_omicron,
                      seed20_39_omicron_ba2, seed40_59_omicron_ba2,
                      seed20_39_omicron_ba4_5, seed40_59_omicron_ba4_5]   # or seed_1_6 if using old way
            );
locdat = popdat[locale];


# %% [markdown]
# ### Some summary statistics

# %%
dead_end = series[locale].cum.dead_total[end] - series[locale].cum.dead_total[390]
dead_beg = series[locale].cum.dead_total[390] - series[locale].cum.dead_total[1]
@show(dead_end, dead_beg)

# %%
stat_cond = cs.stat_cond(series, locale)

# %%
stat_vax = cs.stat_vax(popdat, locale)
# @Select(item, age80_up_pct, age60_79_pct, age40_59_pct, age20_39_pct, age0_19_pct)(stat_vax)

# %%
pretty_table(Tables.matrix(stat_vax, transpose=true)[2:end,:], row_names=collect(columnnames(stat_vax))[2:end], header=stat_vax.item)

# %%
stat_repeat = cs.stat_repeat(popdat, locale)

# %% [markdown]
# ### Serialize model output and modeldef for kicks

# %%
cs.series_to_csv(series, basedir=:none, pathstr="/Users/lewis/Dropbox/Covid Modeling/Outputs")

# %%
cs.popdat_to_csv(popdat, basedir=:none, idstr="eval_vax", pathstr="/Users/lewis/Dropbox/Covid Modeling/Outputs")

# %%
cs.modeldef_to_yaml(ndays, [locale], day1 = Date("2020-01-01", "yyyy-mm-dd"), 
                    dovax=true,
                    basedir=:none, 
                    scheddir = "vaccine_100k",
                    pathstr="/Users/lewis/Dropbox/Covid Modeling/Outputs")

# %% [markdown] tags=[]
# ### Plot results

# %%
cumplot(series, locale, [:unexposed, :infectious, :recovered, :dead, :totvaccinated])

# %% [markdown]
# Note that the orange line labeled Infectious, which shows the current number of infected people, is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: There were some sick people as of the day before. Some more people got sick today. Some people got better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. So net infected is yesterday + new today - recovered today - died today. Newspaper tracking shows the new infections of each day--who got sick today? Tomorrow, if no one new got sick the line would be at zero--even though the people who got sick yesterday are still infected. So, the newspaper line goes up and down faster and higher, because no one is ever subtracted. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line levels off. 

# %%
                  
cs.daily_cases_plot(series, popdat, locale, [:total]; geo=[], thm=:ggplot2)

# %% tags=[]
cumplot(series, locale, [:nil, :mild, :sick, :severe])

# %% tags=[]
newplot(series, locale, [:dead])

# %%
cumplot(series, locale, [:base, :delta, :omicron_ba1, :omicron_ba2, :omicron_ba4_5])


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
stat_cond = cs.stat_cond(series, locale)

# %%
count(locdat.status .== infectious)

# %%
stat_cond.age80_up_pct

# %%
stat_vax = cs.stat_vax(popdat, locale)

# %%
newplot(series, locale, [:Pfizer], days=380:700)


# %%
newplot(series, locale, [:Moderna], days=380:700)

# %%
newplot(series, locale, :JnJ, days=380:700)


# %% [markdown]
# ### Vaccine schedule 

# %%
fieldnames(typeof(model.vaxschedset[:loc38015_old]))

# %%
sf = model.vaxschedset[:loc38015_old].spreadfunc

# %%
model.vaxschedset[:loc38015_old].dayrange

# %%
model.vaxschedset[:loc38015_old].targetpct

# %%
sum(sf.(350:700))

# %%
plot(sf.(350:700))

# %%
sum(sf.(350:2:700))  
# if you don't query all the days, you can get much less than the full sum--but we do query all days
# if you don't use all the allocation on a "big" day, you can't get that quantity back
# so what's the solution:
      # only impose a limit on the early days
      # then don't enforce any limit later
      # this can be easier than using the interpolation generator funcion
      # the called function should make the date range more obvious
      # the only consequence is running out early, but that is what would really happen
# is it better to have a limit on people or doses?
      # a limit on people is modeling the behavior of requesting the virus
      # a limit on doses is the real world of limited supply of vaccine
      # we are not going to track the number of people turned away, because 1) it doesn't matter; 2) the number generated
          # by the simulation has no relation to a similar problem in the real world
# should the limit be calculated on the no. of doses currently available or the original starting supply?
      # pct of what's left means we never use it up
      # pct of original supply means we never tap previous days unused doses
      # might be simpler to use doses / day such that all the vaccine is used up
      # we could also have an explict "request" schedule, especially for the early weeks

# %%
struct sched
    pattern::Vector{Float64}
    startday::Int
end

# %%
seq = [0.02, .05, .08, .15, .2, .24, .05, .05, .03, .03, .03, .03, .02, .02]

# %%
s1 = sched(seq, 50)

# %%
startqty = 1000

function dispense(startq, sch, day)
    if day < sch.startday
        return 0
    end
    weeknum = fld(day - sch.startday, 7) + 1
    #  portion = weeknum <= length(sch.pattern) ? sch.pattern[weeknum] : sch.pattern[end]
    weekportion = sch.pattern[max(weeknum, length(sch.pattern))]
    give = round(Int, (startq * weekportion) / 7, RoundDown)
end

# %%
dispense(5_000, s1, 52)

# %%
model.vaxschedset[:loc38015_old].vaxesincluded

# %%
model.vaxschedset[:loc38015_old].vaxesincluded[:Moderna]

# %%
newplot(series, locale, :unexposed)




# %% [markdown]
# ## How many people got sick multiple times?
#
#

# %%
sort(countmap(length.(locdat.variant)))

# %% [markdown]
# ### Breakout Infections

# %%
@Select(vaxday, sickday, agegrp)(locdat[((locdat.status .== infectious) .| (locdat.status .== recovered)) .& 
                            (locdat.vaxstatus .!= :none) .& (last.(locdat.sickday) .> first.(locdat.vaxday))])

# %% [markdown]
# ### How many vaccinated people died?

# %%
@Select(agegrp, vaxday, sickday, deadday)(locdat[(locdat.status .== dead) .& (locdat.vaxstatus .!= :none)])

# %% [markdown]
# ## Test a social distancing case

# %%
sd1 = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    

# %%
sd1_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, include_ages=[])

# %%
popdat, series = runsim(model;
            dovax=true, vaxscheds=:all, showr0=false, silent=true, 
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, seed20_39_omicron, seed40_59_omicron,
                      seed20_39_omicron_ba2, seed40_59_omicron_ba2, seed20_39_omicron_ba4_5, 
                      seed40_59_omicron_ba4_5, sd1, sd1_end]);
locdat = popdat[locale];

# %%
cumplot(series, locale, [:unexposed, :infectious, :recovered, :dead, :totvaccinated])

# %%
cumplot(series, locale,[:infectious, :dead])

# %%
cumplot(series, locale,[:base, :delta, :omicron_ba1, :omicron_ba2, :omicron_ba4_5])

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
                        dovax=true, vaxscheds=:all,
                        runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, 
                                  seed20_39_omicron, seed40_59_omicron,
                                  seed20_39_omicron_ba2, seed40_59_omicron_ba2, seed20_39_omicron_ba4_5, seed40_59_omicron_ba4_5,
                                  sdolder, sdolder_end]);
locdat = popdat[locale];


# %%
olderdat =popdat[locale]

sd = findall(olderdat.sdcase .!= :none)

@Select(agegrp, sdcase)(olderdat[sd])

# %%
cumplot(series, locale, [:unexposed, :infectious, :recovered, :dead, :totvaccinated])

# %%
cumplot(series, locale, [:infectious, :dead])

# %% [markdown]
# ## Social Distancing starts with everyone and then the younger folks party

# %%
sdyoung_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age0_19, age20_39])    

# %%
popdat, series = runsim(model,
                        dovax=true, vaxscheds=:all,
                        runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, 
                                  seed20_39_omicron, seed40_59_omicron, seed20_39_omicron_ba2, 
                                  seed40_59_omicron_ba2, seed20_39_omicron_ba4_5, seed40_59_omicron_ba4_5,
                                  sd1, sdyoung_end]);
locdat = popdat[locale];

# %%
mixdat = popdat[locale]

sd = findall(mixdat.sdcase .!= :none)
young = findall((mixdat.agegrp .== age0_19) .| (mixdat.agegrp .== age20_39))
sd_young_idx = intersect(sd, young)
old = findall((mixdat.agegrp .== age40_59) .| (mixdat.agegrp .== age60_79) .| (mixdat.agegrp .== age80_up));

# %%
youngtab = Table(mixdat[young])
count(youngtab.sdcase .== :none)

# %%
oldtab = Table(mixdat[old])
count(oldtab.sdcase .!= :none)

# %%
typeof(mixdat.agegrp)

# %%
mixdat = popdat[locale]

# %%
cumplot(series, locale, [:unexposed, :infectious, :recovered, :dead, :totvaccinated])

# %%
@Select(status, agegrp, cond, sdcomply)(locdat)
