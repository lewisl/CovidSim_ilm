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
#     display_name: Julia 1.7.2
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
using Distributions
using YAML
using PrettyPrint
using LinearAlgebra
using Dates

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# # Test setup and population matrix

# %% tags=[]
ndays = 180
locale = 38015
model180 = buildsim(ndays, locale;  
    day1 = Date("2020-01-01", "yyyy-mm-dd"),
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
println(typeof(age80tree))
breakday_idx = 5
println("for breakday $breakday_idx")
age80tree
# println("fieldnames: ", fieldnames(typeof(age80tree[breakday_idx])))
# println("duration ",age80tree[breakday_idx].duration)
# println("transition \n", age80tree[breakday_idx].transition)

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
ndays = 850
locale = 38015
model = buildsim(ndays, locale;  
    day1 = Date("2020-01-01", "yyyy-mm-dd"),
    dovax = true,
    paramdir = "../parameters",
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
            dovax=true,
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_delta, seed40_59_delta, seed20_39_omicron, seed40_59_omicron,
                      seed20_39_omicron_ba2, seed40_59_omicron_ba2]   # or seed_1_6 if using old way
            );
locdat = popdat[locale];


# %% [markdown]
# ### Plot results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orange line labeled Infectious, which shows the current number of infected people, is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: There were some sick people as of the day before. Some more people got sick today. Some people got better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. So net infected is yesterday + new today - recovered today - died today. Newspaper tracking shows the new infections of each day--who got sick today? Tomorrow, if no one new got sick the line would be at zero--even though the people who got sick yesterday aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. 

# %% jupyter={"outputs_hidden": true} tags=[]
cumplot(series, locale, [:nil, :mild, :sick, :severe, :totinfected])

# %% jupyter={"outputs_hidden": true} tags=[]
newplot(series, locale, [:dead])

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
