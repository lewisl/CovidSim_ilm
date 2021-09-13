# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %%
using CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML
using PrettyPrint
using Plots

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% [markdown]
# # Test setup and population matrix

# %% tags=[]
# set locale and number of days
locale = 38015
ndays = 180

# %% tags=[]
alldat = setup(ndays, [locale]; 
    dovax=true,
    paramdir="../parameters", 
    geofilename="../data/geo2data.csv",
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml"
    );

# %%
keys(alldat)

# %%
alldat.dat

# %%
locdat = alldat.dat["popdat"][locale]

# %%
locdat.vaxrcvd

# %%
ages = alldat.dat["agegrp_idx"][locale]

# %%
columnnames(locdat)

# %%
countmap(locdat.agegrp)

# %%
countmap(locdat.status)  # everyone begins as unexposed

# %% tags=[]
geodf = alldat.geo   # the date for all locales has been read into a dataframe

# %%
density_factor = geodf[geodf[!, :fips] .== locale, :density_factor][]

# %%
alldat.spreadset[:default]  # the spread parameters are loaded as a dict of float arrays

# %%
alldat.social

# %%
fieldnames(typeof(alldat.social))

# %%
contactfactors = alldat.social.contactfactors

# %%
typeof(contactfactors)

# %%
contactfactors[age80_up]

# %%
touchfactors =  alldat.social.touchfactors

# %%
touchfactors[age40_59]

# %%
limdict = CovidSim_ilm.limdict
limdict(touchfactors, <)  # recursive minimum

# %%
alldat.vaxset

# %%
# is shifter working?
shifter(touchfactors, (.18, .3)...)[age40_59]

# %%
dectree = alldat.transitionset[:default] # the decision trees for all age groups are loaded

# %%
transarr = alldat.transitionset[:default]

# %%
transarr[age0_19]

# %%
transarr[age0_19][1]

# %%
transarr[age0_19][1][:transition]

# %%
tr_age0_19 = transarr[age0_19]
T = typeof(tr_age0_19[1][:transition])
supertype(T)

# %%
trvec = zeros(6)
p_cond = mild
p_sickday = 9
trvec = CovidSim_ilm.has(tr_age0_19, p_sickday, p_cond, trvec)
println(typeof(trvec))
trvec



# %%
typeof(dectree)

# %% [markdown]
# Dict{Int64, OrderedCollections.OrderedDict{Int, Dict{String, Vector{T} where T}

# %% tags=[]
dectree[age80_up]

# %%
typeof(dectree[age80_up][5][:transition])

# %% [markdown]
# # Load vaccine parameters

# %%
# mapping for all vaccines to parameters for each vaccine
vaccines = YAML.load_file("../parameters/vaccines.yml"; dicttype=Dict{Symbol,Any})

# %%
vaccines[:Pfizer]

# %%
vaxkeys = keys(vaccines)

# %%
vaccines[:Pfizer]

# %% [markdown]
# ### Other vaccine parameters

# %%
vax = :Moderna
v = YAML.load_file(joinpath("../parameters","vaccine_parameters",vaccines[vax][:directory_name],
        vaccines[vax][:infect_fname]), dicttype=Dict{Symbol, Any})

# %%
vaxset = CovidSim_ilm.build_vaxset("vaccines.yml", paramdir="../parameters")

# %% [markdown]
# ### Vaccination Schedule

# %%
vxschedset = CovidSim_ilm.build_vaxschedset()


# %%
println(typeof(vxschedset))
pprint(vxschedset)

# %%
asched = first(keys(vxschedset))

# %%
vxschedset[asched].dayrange

# %%
sum([vxschedset[asched].pctfunc(i) for i in vxschedset[asched].dayrange])

# %%
plot(vxschedset[asched].dayrange, [vxschedset[asched].pctfunc(i) for i in vxschedset[asched].dayrange],size=(600,300))

# %% [markdown] tags=[]
# # Give some shots

# %%
peeps = 30:40  
cnt = length(peeps)
day_ctr[:day] = 625  # must be in the range of the schedule
vaccinate!(locdat, vxschedset, peeps, vaxset)

# %%
@Select(vaxstatus, vaxday, vaxrcvd)(locdat)[peeps]

# %% [markdown]
# #### Reset the give shots test

# %%
locdat.vaxstatus[peeps] = fill(:none, cnt)
locdat.vaxday[peeps] = fill([0], cnt)
locdat.vaxrcvd[peeps] = fill([:none], cnt)

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :default, agegrps)

# %% [markdown]
# # Run a simulation

# %%
result_dict, series = run_a_sim(ndays, locale; 
    dovax=false, 
    dovariant=false,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml",
    showr0=false, 
    silent=true, 
    runcases=[seed_1_6]);

# %%
keys(result_dict)

# %%
keys(result_dict[:dat])

# %%
popdat = result_dict[:dat]["popdat"][locale]

# %%
countmap(popdat.cond)

# %%
countmap(popdat.status)

# %%
virus_outcome(series, locale, base=:pop)

# %%
series[locale][:cum]

# %% [markdown]
# # Plot results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.

# %% [markdown]
# ## Run simulation with full vaccination schedule

# %%
ndays = 720
result_dict, series = run_a_sim(ndays, locale; 
    dovax=true, 
    dovariant=false,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml",
    showr0=false, 
    silent=true, 
    runcases=[seed_1_6]);

# %%
locdat = result_dict[:dat]["popdat"][locale]
vaxcols = @Select(vaxrcvd, vaxstatus)(locdat)

# %%
pfizer2shots = sum(map(v -> (first(v) == :Pfizer) & (length(v) == 2) ? 2 
            : 0, vaxcols.vaxrcvd))
pfizer1shots = sum(map(v -> (first(v) == :Pfizer) & (length(v) == 1) ? 1 
            : 0, vaxcols.vaxrcvd))

peoplepfizer2 = sum(map(v -> (first(v) == :Pfizer) & (length(v) == 2) ? 1 
            : 0, vaxcols.vaxrcvd))
peoplepfizer1 = sum(map(v -> (first(v) == :Pfizer) & (length(v) == 1) ? 1 
            : 0, vaxcols.vaxrcvd))

moderna2shots = sum(map(v -> (first(v) == :Moderna) & (length(v) == 2) ? 2 
            : 0, vaxcols.vaxrcvd))
moderna1shots = sum(map(v -> (first(v) == :Moderna) & (length(v) == 1) ? 1 
            : 0, vaxcols.vaxrcvd))

peoplemoderna2 = sum(map(v -> (first(v) == :Moderna) & (length(v) == 2) ? 1 
            : 0, vaxcols.vaxrcvd))
peoplemoderna1 = sum(map(v -> (first(v) == :Moderna) & (length(v) == 1) ? 1 
            : 0, vaxcols.vaxrcvd))

peoplejnj = jnjshots = sum(map(v -> (first(v) == :JnJ) & (length(v) == 1) ? 1 
            : 0, vaxcols.vaxrcvd))

peopleatleast1shot = (peoplepfizer2 + peoplepfizer1 + peoplemoderna2 + peoplemoderna1 +
                    peoplejnj)

peoplefull = count(vaxcols.vaxstatus .== :full) 

@show pfizer2shots, pfizer1shots
@show peoplepfizer2, peoplepfizer1

@show moderna2shots, moderna1shots
@show peoplemoderna2, peoplemoderna1

@show jnjshots

@show peopleatleast1shot, peoplefull

# %%
count(vaxcols.vaxstatus .== :full)

# %%
count(vaxcols.vaxstatus .== :first)

# %% [markdown]
# ## Test a social distancing case

# %%
sd1 = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    

# %% tags=[]
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
outdat = result_dict[:dat]["popdat"][locale]
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
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, 
    runcases=[seed_1_6, sdolder, sdolder_end]);


# %%
olderdat = result_dict[:dat]["popdat"][locale]

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
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, 
    runcases=[seed_1_6, sd1, sdyoung_end]);

# %%
mixdat = result_dict[:dat]["popdat"][locale]

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
mixdat = result_dict[:dat]["popdat"][locale]

# %%
cumplot(series, locale, [:infectious, :dead])

# %%
@Select(status, agegrp, cond, sdcomply)(locdat)

# %% [markdown]
# alldat

# %% tags=[]
alldat

# %% tags=[]
ages = alldat.dat["agegrp_idx"][locale]

# %%
include_ages = [age0_19, age20_39]

# %%
union((ages[i] for i in include_ages)...)

# %%
locdat.sdcomply[collect(1:5:95000)] .= :test

# %%
incase_idx = findall(locdat.sdcomply .== :test)

# %%
byage_idx = intersect(incase_idx, union((ages[i] for i in include_ages)...))

# %%
@btime locdat.status;

# %%
statuscol = locdat.status
@btime statuscol;

# %%
nt = (one=1, two=2, three=3)

# %%
@btime nt.one;

# %%
thisone = nt.one
@btime thisone;

# %%
