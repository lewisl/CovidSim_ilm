# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% jupyter={"outputs_hidden": true}
pwd()

# %% jupyter={"outputs_hidden": true}
using CovidSim_ilm

# %% jupyter={"outputs_hidden": true}
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML
using PrettyPrint
using Plots
# Plots.pyrcparams["backend"]="Qt5Agg"

# %% [markdown]
# # Test setup and population matrix

# %% tags=[] jupyter={"outputs_hidden": true}
# set locale and number of days
locale = 38015
ndays = 180

# %% jupyter={"outputs_hidden": true}
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# %% tags=[] jupyter={"outputs_hidden": true}
alldat = setup(ndays, [locale]; 
    dovax=true,
    paramdir="../parameters", 
    geofilename="../data/geo2data.csv",
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml"
    );

# %% jupyter={"outputs_hidden": true}
keys(alldat)

# %% jupyter={"outputs_hidden": true}
alldat.dat

# %% jupyter={"outputs_hidden": true}
locdat = alldat.dat["popdat"][locale]

# %% jupyter={"outputs_hidden": true}
locdat.vaxrcvd

# %% jupyter={"outputs_hidden": true}
ages = alldat.dat["agegrp_idx"][locale]

# %% jupyter={"outputs_hidden": true}
columnnames(locdat)

# %% jupyter={"outputs_hidden": true}
countmap(locdat.agegrp)

# %% jupyter={"outputs_hidden": true}
countmap(locdat.status)  # everyone begins as unexposed

# %% tags=[] jupyter={"outputs_hidden": true}
geodf = alldat.geo   # the date for all locales has been read into a dataframe

# %% jupyter={"outputs_hidden": true}
density_factor = geodf[geodf[!, :fips] .== locale, :density_factor][]

# %% jupyter={"outputs_hidden": true}
pprintln(alldat.spreadset[:base])  # the spread parameters are loaded as a dict of float arrays

# %% jupyter={"outputs_hidden": true}
alldat.social

# %% jupyter={"outputs_hidden": true}
fieldnames(typeof(alldat.social))

# %% jupyter={"outputs_hidden": true}
contactfactors = alldat.social.contactfactors

# %% jupyter={"outputs_hidden": true}
typeof(contactfactors)

# %% jupyter={"outputs_hidden": true}
contactfactors[age80_up]

# %% jupyter={"outputs_hidden": true}
touchfactors =  alldat.social.touchfactors

# %% jupyter={"outputs_hidden": true}
touchfactors[age40_59]

# %% jupyter={"outputs_hidden": true}
limdict = CovidSim_ilm.limdict
limdict(touchfactors, <)  # recursive minimum

# %% jupyter={"outputs_hidden": true}
alldat.vaxset

# %% jupyter={"outputs_hidden": true}
# is shifter working?
shifter(touchfactors, (.18, .3)...)[age40_59]

# %% jupyter={"outputs_hidden": true}
dectree = alldat.transitionset[:base] # the decision trees for all age groups are loaded

# %% jupyter={"outputs_hidden": true}
transarr = alldat.transitionset[:base]

# %% jupyter={"outputs_hidden": true}
transarr[age0_19]

# %% jupyter={"outputs_hidden": true}
transarr[age0_19][1]

# %% jupyter={"outputs_hidden": true}
transarr[age0_19][1][:transition]

# %% jupyter={"outputs_hidden": true}
tr_age0_19 = transarr[age0_19]
T = typeof(tr_age0_19[1][:transition])
supertype(T)

# %% jupyter={"outputs_hidden": true}
trvec = zeros(6)
p_cond = mild
p_sickday = 9
trvec = CovidSim_ilm.has(tr_age0_19, p_sickday, p_cond, trvec)
println(typeof(trvec))
trvec



# %% jupyter={"outputs_hidden": true}
typeof(dectree)

# %% [markdown]
# Dict{Int64, OrderedCollections.OrderedDict{Int, Dict{String, Vector{T} where T}

# %% tags=[] jupyter={"outputs_hidden": true}
dectree[age80_up]

# %% jupyter={"outputs_hidden": true}
typeof(dectree[age80_up][5][:transition])

# %% [markdown]
# ## Let's look at the decay rate of immunity gained by recovering

# %% jupyter={"outputs_hidden": true}
CovidSim_ilm.risk(alldat.spreadset, alldat.vaxset, locdat, 20, 200)
spreadset, vaxset, locdat, spreader, target

# %% [markdown]
# # Load vaccine parameters

# %% jupyter={"outputs_hidden": true}
# mapping for all vaccines to parameters for each vaccine
vaccines = YAML.load_file("../parameters/vaccines.yml"; dicttype=Dict{Symbol,Any})

# %% jupyter={"outputs_hidden": true}
vaccines[:Pfizer]

# %% jupyter={"outputs_hidden": true}
vaxkeys = keys(vaccines)

# %% jupyter={"outputs_hidden": true}
vaccines[:Pfizer]

# %% [markdown]
# ### Other vaccine parameters

# %% jupyter={"outputs_hidden": true}
vax = :Moderna
v = YAML.load_file(joinpath("../parameters","vaccine_parameters",vaccines[vax][:directory_name],
        vaccines[vax][:infect_fname]), dicttype=Dict{Symbol, Any})

# %% jupyter={"outputs_hidden": true}
vaxset = CovidSim_ilm.build_vaxset("vaccines.yml", paramdir="../parameters")

# %% [markdown]
# ### Vaccination Schedule

# %% jupyter={"outputs_hidden": true}
vxschedset = CovidSim_ilm.build_vaxschedset()


# %% jupyter={"outputs_hidden": true}
println(typeof(vxschedset))
pprint(vxschedset)

# %% jupyter={"outputs_hidden": true}
asched = first(keys(vxschedset))

# %% jupyter={"outputs_hidden": true}
vxschedset[asched].dayrange

# %% jupyter={"outputs_hidden": true}
sum([vxschedset[asched].pctfunc(i) for i in vxschedset[asched].dayrange])

# %% jupyter={"outputs_hidden": true}
plot(vxschedset[asched].dayrange, [vxschedset[asched].pctfunc(i) for i in vxschedset[asched].dayrange],size=(600,300))

# %% [markdown] tags=[]
# # Give some shots

# %% jupyter={"outputs_hidden": true}
peeps = 301:400
cnt = length(peeps)
day_ctr[:day] = 500  # must be in the range of the schedule
vaccinate!(locdat, vxschedset, peeps, vaxset)

# %% jupyter={"outputs_hidden": true}
any(last.(locdat.vaxrcvd) .!= :none)

# %% [markdown]
# #### Reset the give shots test

# %% jupyter={"outputs_hidden": true}
locdat.vaxstatus[peeps] = fill(:none, cnt)
locdat.vaxday[peeps] = fill([0], cnt)
locdat.vaxrcvd[peeps] = fill([:none], cnt)

# %% [markdown]
# # Create a seed case

# %% jupyter={"outputs_hidden": true}
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :base, agegrps)

# %% [markdown]
# # Run a simulation

# %% tags=[]
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

# %% tags=[]
keys(result_dict)

# %% tags=[]
rf(x) = CovidSim_ilm.sigmoid(shifter(clamp(x, 0.0, 2.0), 0.0, 2.0, -4.5, 4.5))  

# %% tags=[]
rf(2.0)

# %% tags=[]
keys(result_dict[:dat])

# %% tags=[]
popdat = result_dict[:dat]["popdat"][locale]

# %% tags=[]
countmap(popdat.cond)

# %% tags=[]
countmap(popdat.status)

# %% tags=[]
virus_outcome(series, locale, base=:pop)

# %% tags=[]
series[locale][:cum]

# %% [markdown]
# # Plot results

# %% tags=[]
cumplot(series, locale)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.

# %% [markdown]
# ## Run simulation with full vaccination schedule

# %% jupyter={"outputs_hidden": true}
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

# %% jupyter={"outputs_hidden": true}
locdat = result_dict[:dat]["popdat"][locale]
vaxcols = @Select(vaxrcvd, vaxstatus)(locdat)

# %% jupyter={"outputs_hidden": true}
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

# %% jupyter={"outputs_hidden": true}
count(vaxcols.vaxstatus .== :full)

# %% jupyter={"outputs_hidden": true}
count(vaxcols.vaxstatus .== :first)

# %% [markdown]
# ## Test a social distancing case

# %% jupyter={"outputs_hidden": true}
sd1 = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    

# %% tags=[] jupyter={"outputs_hidden": true}
sd1_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, include_ages=[])

# %% jupyter={"outputs_hidden": true}
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, runcases=[seed_1_6, sd1, sd1_end]);

# %% jupyter={"outputs_hidden": true}
virus_outcome(series, locale, base=:pop)

# %% jupyter={"outputs_hidden": true}
cumplot(series, locale)

# %% jupyter={"outputs_hidden": true}
cumplot(series, locale,[:infectious, :dead])

# %% jupyter={"outputs_hidden": true}
outdat = result_dict[:dat]["popdat"][locale]
all(outdat.sdcomply .== :none)

# %% [markdown]
# ## Social distancing only among those age40_59, age60_79, age80_plus

# %% jupyter={"outputs_hidden": true}
sdolder = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# %% jupyter={"outputs_hidden": true}
sdolder_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# %% jupyter={"outputs_hidden": true}
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, 
    runcases=[seed_1_6, sdolder, sdolder_end]);


# %% jupyter={"outputs_hidden": true}
olderdat = result_dict[:dat]["popdat"][locale]

sd = findall(olderdat.sdcomply .!= :none)

@Select(agegrp, sdcomply)(olderdat[sd])

# %% jupyter={"outputs_hidden": true}
cumplot(series, locale)

# %% jupyter={"outputs_hidden": true}
cumplot(series, locale, [:infectious, :dead])

# %% [markdown]
# ## Social Distancing starts with everyone and then the younger folks party

# %% jupyter={"outputs_hidden": true}
sdyoung_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age0_19, age20_39])    

# %% jupyter={"outputs_hidden": true}
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, 
    runcases=[seed_1_6, sd1, sdyoung_end]);

# %% jupyter={"outputs_hidden": true}
mixdat = result_dict[:dat]["popdat"][locale]

sd = findall(mixdat.sdcomply .!= :none)
young = findall((mixdat.agegrp .== age0_19) .| (mixdat.agegrp .== age20_39))
sd_young_idx = intersect(sd, young)
old = findall((mixdat.agegrp .== age40_59) .| (mixdat.agegrp .== age60_79) .| (mixdat.agegrp .== age80_up));

# %% jupyter={"outputs_hidden": true}
youngtab = Table(mixdat[young])
count(youngtab.sdcomply .== :none)

# %% jupyter={"outputs_hidden": true}
oldtab = Table(mixdat[old])
count(oldtab.sdcomply .!= :none)

# %% jupyter={"outputs_hidden": true}
typeof(mixdat.agegrp)

# %% jupyter={"outputs_hidden": true}
mixdat = result_dict[:dat]["popdat"][locale]

# %% jupyter={"outputs_hidden": true}
cumplot(series, locale, [:infectious, :dead])

# %% jupyter={"outputs_hidden": true}
@Select(status, agegrp, cond, sdcomply)(locdat)

# %% [markdown]
# alldat

# %% tags=[] jupyter={"outputs_hidden": true}
alldat

# %% tags=[] jupyter={"outputs_hidden": true}
ages = alldat.dat["agegrp_idx"][locale]

# %% jupyter={"outputs_hidden": true}
include_ages = [age0_19, age20_39]

# %% jupyter={"outputs_hidden": true}
union((ages[i] for i in include_ages)...)

# %% jupyter={"outputs_hidden": true}
locdat.sdcomply[collect(1:5:95000)] .= :test

# %% jupyter={"outputs_hidden": true}
incase_idx = findall(locdat.sdcomply .== :test)

# %% jupyter={"outputs_hidden": true}
byage_idx = intersect(incase_idx, union((ages[i] for i in include_ages)...))

# %% jupyter={"outputs_hidden": true}
@btime locdat.status;

# %% jupyter={"outputs_hidden": true}
statuscol = locdat.status
@btime statuscol;

# %% jupyter={"outputs_hidden": true}
nt = (one=1, two=2, three=3)

# %% jupyter={"outputs_hidden": true}
@btime nt.one;

# %% jupyter={"outputs_hidden": true}
thisone = nt.one
@btime thisone;

# %% jupyter={"outputs_hidden": true}
