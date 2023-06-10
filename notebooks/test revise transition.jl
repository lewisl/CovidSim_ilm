# ---
# jupyter:
#   jupytext:
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
export @ballocated, @belapsed, @benchmark, @benchmarkable, @benchmarkset, @bprofile, @btime, @case, @tagged, BenchmarkGroup, BenchmarkTools, addgroup!, allocs, gctime, improvements, invariants, isimprovement, isinvariant, isregression, judge, leaves, loadparams!, mean, median, memory, params, ratio, regressions, rmskew, rmskew!, trim, tune!, warmup
using Distributions
using YAML
using PrettyPrint

# %%
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/src"))

# %% [markdown]
# # Test setup and population matrix

# %% tags=[]
# set locale and number of days
locale = 38015
ndays = 180

# %% tags=[]
alldat = setup(ndays, [locale]; paramdir="../parameters", 
    geofilename="../data/geo2data.csv",
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml"
    );

# %%
alldat.transitionset # the decision transition matrices for all age groups are loaded

# %%
dectree = alldat.transitionset[:base]

# %%
dectree[age80_up]

# %%
dectree[age80_up][5]

# %%
trarr = copy(dectree[age80_up][5][:transition])

# %%
axes(trarr, 1)

# %%
function adjust_transition!(trarr, delta= 0.05)
    r, c = size(trarr)
    half = round(Int, c/2)
    adjvec = reshape(zeros(c), 1, c)
    prev_val = 0
    for i = 1:c
        if i <= half
            val = 1.0 + float(i) * delta
            prev_val = val
        else
            val = prev_val - float(i) * delta
        end
        adjvec[1, i] = val
    end

    for i in axes(trarr, 1)
        rowsum = sum(trarr[i,:])
        if rowsum != 0.0
            for c in axes(trarr,2)
                trarr[i, c] *= adjvec[c]
            end
            nf = 1.0 / sum(trarr[i,:])
            trarr[i,:] .*= nf
        end
    end
end

# %%
function adjust_transition2!(trarr; vaxtreffect= 0.9, vaxdiscount=0.5)
    for i in axes(trarr, 1)
        rowsum = sum(trarr[i,:])
        if rowsum != 0.0
            for c in (sick, severe, dead)
                trarr[i, CovidSim_ilm.map2outcome(c)] *= (1.0 - vaxtreffect)
            end
            nf = 1.0 / sum(trarr[i,:])
            trarr[i,:] .*= nf
        end
    end
end

# %%
@btime adjust_transition2!($trarr); 

# %%
trarr

# %%
chgvec .* trarr

# %%
sum(chgvec .* trarr, dims=2)

# %%
dectree[age80_up][5][:transition][nil, :]

# %% [markdown]
# # Create a seed case

# %%
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, :nil, :base, agegrps)

# %% [markdown]
# # Run a simulation

# %%
ndays = 180
locale = 38015
model = buildsim(ndays, locale;  
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml",
);

# %%
typeof(model)

# %%
popdat, series = runsim(model;
            runcases=[seed_1_6]
            );

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
virus_outcome(series, locale, base=:pop)

# %%
keys(series[locale])

# %%
println(":cum ", typeof(series[locale][:cum]))
println(":cum ", size(series[locale][:cum]))
println(":new ", size(series[locale][:new]))

# %%
series[locale][:cum]

# %% [markdown]
# # Plotted results

# %%
cumplot(series, locale)

# %% [markdown]
# Note that the orangle line labeled Infectious that shows the number of infected people is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: Some people got sick today. Some people get better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. Newspaper tracking shows the new active infections of each day--who got sick today? The next day, if no one new got sick the line would be at zero--even though the people who got sick aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. This is the least common way to show the data.

# %%
cumplot(series, locale, [:nil, :mild, :sick, :severe])

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
function maketup(keys, values)
    NamedTuple{keys}(values)
end

# %%
mapconds = maketup((:nil, :mild, :sick, :severe), [5,6,7,8])

# %%
typeof(mapconds)

# %%
keys(mapconds)

# %%
eltype(values(mapconds))

# %%
using BenchmarkTools
@btime $mapconds.nil

# %%
function mapvec(maptup, vals)
    [getfield(maptup, x) for x in vals]
end

# %%
typeof(mapvec)

# %%

# %%
vals = [nil, sick]
mapvec(mapconds, Symbol.(vals))

# %%
using Markdown
println(md"*bold*")

# %%
map2series

# %%
getindex(map2series, :unexposed)

# %%
getindex(Dict("a"=>1, "b"=>2, "c"=>3), "a")
