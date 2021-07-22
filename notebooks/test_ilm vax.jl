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
alldat = setup(ndays, [locale])

# %%
alldat.dat

# %%
locdat = alldat.dat["popdat"][locale]

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
infectparams = alldat.infect  # the spread parameters are loaded as a dict of float arrays

# %%
typeof(infectparams)

# %%
socialparams = alldat.social

# %%
fieldnames(typeof(socialparams))

# %%
contactfactors = socialparams.contactfactors

# %%
typeof(contactfactors)

# %%
contactfactors[age80_up]

# %%
touchfactors =  socialparams.touchfactors

# %%
touchfactors[age40_59]

# %%
limdict = CovidSim_ilm.limdict
limdict(touchfactors, <)  # recursive minimum

# %%
# is shifter working?
shifter(touchfactors, (.18, .3)...)[age40_59]

# %%
dectree = alldat["dectree"] # the decision trees for all age groups are loaded

# %%
typeof(dectree)

# %% tags=[]
display_tree(dectree)

# %% [markdown]
# Dict{Int64, OrderedCollections.OrderedDict{Int, Dict{String, Vector{T} where T}

# %% tags=[]
dectree[age80_up]

# %%
typeof(dectree[age80_up][25][sick][:outcomes])

# %%
function get_node(dectree, agegrp, sickday, fromcond)
    dectree[agegrp][sickday][fromcond]
end

# %%
@btime get_node(dectree, age80_up, 25, sick)[:probs]

# %% [markdown]
# # Load vaccine parameters

# %%
# mapping for all vaccines to parameters for each vaccine
vaccines = YAML.load_file("../parameters/vaccines.yml"; dicttype=Dict{Symbol,Any})

# %%
vaxkeys = keys(vaccines)

# %%
vaccines[:Pfizer]

# %% [markdown]
# ### Transition Parameters

# %%
transitionset = Dict()
dectreefilename="../parameters/transition.yml"
for vax in vaxkeys
    transitionset[vax] = setup_dt(joinpath("../parameters", vaccines[vax][:directory_name], 
            vaccines[vax][:transition_fname]))
end
transitionset[:default] = setup_dt(dectreefilename)
transitionset

# %% [markdown]
# ### Other vaccine parameters

# %%
vax = :Moderna
v = YAML.load_file(joinpath("../parameters",vaccines[vax][:directory_name],
        vaccines[vax][:spread_fname]), dicttype=Dict{Symbol, Any})

# %%
spreadset = Dict{Symbol, Union{Infectparams, Vaccineparams}}()
for vax in vaxkeys
    v = YAML.load_file(joinpath("../parameters",vaccines[vax][:directory_name],
        vaccines[vax][:spread_fname]), dicttype=Dict{Symbol, Any})
    v = Vaccineparams(v)
    spreadset[vax] = v
end
spreadset[:default] = Infectparams(; 
    YAML.load_file(joinpath("../parameters","infectparams.yml"), dicttype=Dict{Symbol, Any})...)

# %%
spreadset

# %%
spreadset[:Moderna].recvrisk_reduction

# %%
spreadset[:default]

# %% [markdown]
# ### Vaccination Schedule

# %%
dayrange=60:151
targetpct = 0.75
vaxschedset = Dict{Symbol, Vaxsched}()
for vax in keys(vaccines)
    vaxschedset[vax] = Vaxsched(
        vaccine   = vax,
        dayrange  = dayrange,
        targetpct = targetpct,
        pctperdayfn = makevaxfn(dayrange, spreadset[vax].pattern, targetpct)
        ) 
end


# %%
println(typeof(vaxschedset))
vaxschedset

# %%
vaxschedset[:Pfizer]

# %%
vaxschedset[:Pfizer].pctperdayfn(dayrange.stop)

# %%
sum([vaxschedset[:Pfizer].pctperdayfn(i) for i in dayrange])

# %%
plot(dayrange,[vaxschedset[:Pfizer].pctperdayfn(i) for i in dayrange],size=(600,300))

# %% [markdown] tags=[]
# # Define a vaccine and give some shots

# %%
peeps = 30:40
getashot!(locdat, peeps, 90, :Moderna, vaxkeys)

# %%
locdat.vax[peeps]

# %%
locdat.vaxday[peeps]

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, agegrps)

# %% [markdown]
# # Run a simulation

# %%
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, runcases=[seed_1_6]);

# %%
result_dict

# %%
popdat = result_dict["dat"]["popdat"][locale]

# %%
countmap(popdat.cond)

# %%
countmap(popdat.status)

# %%
virus_outcome(series, locale, base=:pop)

# %%
series[locale][:cum]

# %% [markdown]
# # Plotted results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.

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
outdat = result_dict["dat"]["popdat"][locale]
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
olderdat = result_dict["dat"]["popdat"][locale]

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
mixdat = result_dict["dat"]["popdat"][locale]

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
mixdat = result_dict["dat"]["popdat"][locale]

# %%
cumplot(series, locale, [:infectious, :dead])

# %%
@Select(status, agegrp, cond, sdcomply)(locdat)

# %% [markdown]
# alldat

# %%
alldat

# %%
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
parms = Dict(:one=>1, :two=>2, :v=>[1.2, 2.3])

# %%
Base.@kwdef struct Parms
    one::Int
    two::Int
    v::Vector{Float64}
end

# %%
p1 = Parms(one=1, two=2, v=[1.0, 2.0])

# %%
p1.v

# %%
p2=Parms(;parms...)

# %%
p2.v

# %%
arr = [:alpha, :beta, :gamma, :delta, nothing]

# %%
@btime findfirst(isequal(:beta), arr)

# %%
@btime findfirst(arr .== :beta)

# %%
@btime findfirst(x->x==:beta, arr)

# %%
@btime indexin([:beta], arr)[]

# %%
function findit(item, arr)
    i = 0
    found = false
    for it in arr
        i+=1
        if item == it
            found = true
            break
        end
    end
    return i
end

# %%
@btime begin; i = findit(:delta, arr); arr[i]; end

# %%
m = Dict("alpha"=>0.98, "beta"=>0.94, "gamma"=>0.9, "delta"=>0.84)

# %%
m = Dict(zip(Symbol.(keys(m)), values(m)))

# %%
@btime m[:beta]

# %%
nt = (alpha=1, beta=2, gamma=3, delta=4)

# %%
@btime nt.beta

# %%
@btime getindex(nt, :gamma)

# %%
