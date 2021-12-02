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

# %%
push!(LOAD_PATH, joinpath(homedir(), "Dropbox/Covid Modeling/Covid-ILM/source"))
using CovidSim_ilm

# %%
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML
using PrettyPrint
using Plots
using SplitApplyCombine

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
alldat.spreadset[:base]  # the spread parameters are loaded as a dict of float arrays

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
pprintln(alldat.vaxset[:Pfizer][:params])

# %%
# is shifter working?
shifter(touchfactors, (.18, .3)...)[age40_59]

# %% [markdown]
# ### Testing/checking decision trees (new type as transition matrices)

# %%
dectree = alldat.transitionset[:base] # the decision trees for all age groups are loaded

# %%
transarr = alldat.transitionset[:base]

# %% tags=[]
CovidSim_ilm.display_tree_array(transarr)

# %%
transarr[age0_19]

# %% tags=[]
transarr[age0_19][1][:transition]

# %%
sequences = CovidSim_ilm.getseqs_array(transarr[age0_19]);

# %%
CovidSim_ilm.verifyprobs(sequences)

# %%
CovidSim_ilm.sanitycheck_array(transarr)

# %%
next_transition = [0.9 0.0 0.0 0.1 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.95 0.05 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]

# %%
next_steps = [ i for i in eachindex(next_transition[:,1]) if any(next_transition[i,:] .!= 0.0) ] 

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
transitions = dectree[age80_up][5][:transition]

# %%
breakdays = [dectree[age80_up][i][:sickday] for i in sort(collect(keys(dectree[age80_up])))]

# %%
next_steps = [i for i in eachindex(transitions[:,1]) if any(transitions[i,:] .!= 0.0)]  # if any(transitions[i,:] .!= 0.0) 

# %%
breakday = 5
findfirst(isequal(breakday), breakdays)

# %%
findall(transitions[7,:] .!= 0.0)

# %%
collect(eachindex(transitions[5,:]))

# %%
typeof(dectree[age80_up][5][:transition])

# %%
typeof(dectree[age80_up][5][:transition]) <: AbstractArray

# %%
dectree[age80_up][1][:transition]

# %%
typeof(dectree[age80_up][5][:transition])

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
# ### Other vaccine parameters

# %%
vax = :Moderna
v = YAML.load_file(joinpath("../parameters","vaccine_parameters",vaccines[vax][:directory_name],
        vaccines[vax][:infect_fname]), dicttype=Dict{Symbol, Any})
pprintln(v)

# %%
vaxset = CovidSim_ilm.build_vaxset("vaccines.yml", paramdir="../parameters")

# %% [markdown]
# ### Vaccination Schedule

# %%
vxschedset = CovidSim_ilm.build_vaxschedset()


# %%
println(typeof(vxschedset))
pprintln(vxschedset)

# %%
asched = first(keys(vxschedset))

# %%
targetpct=0.65
pattern=  [0.0, 
          0.02, 
          0.05, 
          0.1, 
          0.15, 
          0.19, 
          0.21, 
          0.16, 
          0.08, 
          0.03, 
          0.01, 
          0.0]

# %%
sum(pattern)

# %%
vxschedset[asched].dayrange

# %%
sum([vxschedset[asched].pctfunc(i) for i in vxschedset[asched].dayrange])

# %%
plot(vxschedset[asched].dayrange, [vxschedset[asched].pctfunc(i) for i in vxschedset[asched].dayrange],size=(600,300))

# %% [markdown] tags=[]
# # Give some shots

# %%
peeps = 1000:1999
cnt = length(peeps)
day_ctr[:day] = 555  # must be in the range of the schedule
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
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :base, agegrps)

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
countmap(popdat.status)

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

# %% [markdown]
# ## Creating filter and seed

# %%
colfetch = :agegrp
op = ==

# %%
filtdict = Dict(:agegrp => [age0_19, age80_up],
                :status => [unexposed],)
filtcmp = Dict(:agegrp => ==, 
                :variant => ==, 
                :cond => ==, 
                :status => ==,
                :sd_comply => ==)

# %%
@btime andor = .&(  
              .|( .==(locdat.agegrp, age0_19), .==(locdat.agegrp, age80_up) ), 
              .|( .==(locdat.status, unexposed), .==(locdat.status, infectious) )
           );

# %%
andor = .&(  
              .|( .==(locdat.agegrp, age0_19), .==(locdat.agegrp, age80_up) ), 
              .|( .==(locdat.status, unexposed), .==(locdat.status, infectious) )
           );
count(andor)

# %%
count(reduce(.|,map(x -> .==(locdat.agegrp, x), filtdict[:agegrp])))

# %%
infilt = [ :agegrp => [age80_up, age0_19], :status => [unexposed, infectious]]

# %%
colage = getproperty(locdat, :agegrp)
colstatus = getproperty(locdat, :status)
@btime ((colage .== age0_19) .| (colage .== age80_up)) .& ((colstatus .== unexposed) .| (colstatus .== infectious));

# %%
@btime reduce(.&, map(y -> reduce(.|, Iterators.map(x -> getproperty($locdat, y.first) .== x, y.second)), $infilt));

# %%
e = :(10 + 10)

# %%
eval(e)

# %%
@eval $e

# %%
holder = [BitVector(rand(Bool, 5)), BitVector(rand(Bool, 5))]

# %%
typeof(locdat.agegrp .== age0_19)

# %% tags=[]
function dofilt_mapreduce(dat, filts)
    # println(filts)
    res = trues(size(dat,1));
    for filt in filts
        comparevals = filt[2]
        col = filt[1]
        res[:] .&= mapreduce(y -> .==(getproperty(dat, col), y), .|, comparevals)
    end
    return res
end

# %% tags=[]
function dofilt_loop(dat, filts)
    # println(filts)
    res = trues(size(dat,1));
    for filt in filts
        comparevals = filt[2]
        col = filt[1]
        for cval in comparevals
            getproperty(dat, col) .== cval
        end
        res[:] .&= mapreduce(y -> .==(getproperty(dat, col), y), .|, comparevals)
    end
    return res
end

# %% tags=[]
finalfilt = dofilt_loop(locdat, infilt);
count(finalfilt)

# %%
finalfilt = dofilt_mapreduce(locdat, infilt);
count(finalfilt)

# %%
@btime dofilt_mapreduce(locdat, infilt);

# %%
@btime dofilt_loop(locdat, infilt);

# %%
@btime foo = ((locdat.agegrp .== age0_19) .| (locdat.agegrp .== age80_up)) .& ((locdat.status .== unexposed) .| (locdat.status .== infectious));

# %%
@btime locdat.status;

# %%
@btime getproperty(locdat, :status);

# %% [markdown]
# ## Using Split Apply Combine

# %%
?group

# %%
selcols = Table(agegrp = locdat.agegrp, status = locdat.status)

# %%
by2colidx = groupinds(selcols)

# %%
locdat[by2colidx[(agegrp=age0_19, status=recovered)]]

# %%
typeof(locdat)

# %%
Popdat = typeof(locdat)

# %%
function getter(dat::Popdat, prop::Symbol, p)
    getproperty(dat, prop)[p]
end

# %%
@btime getter(locdat, :status, 1);

# %%
@btime locdat.status[1];

# %%
col_status = locdat.status;

# %%
@btime col_status[1];

# %%
const a,b,c,d = 1,2,3,4

# %%
foobar = [[1,2,3,4],[10,11,12,13]]

# %%
@enum Colnums begin
    c1 = 1
    c2
    c3
    c4
end

# %%
@btime foobar[Int(c1)][2];

# %%
foobar1 = hcat(foobar...)

# %%
@btime foobar1[1,1];

# %%
@btime foobar[1][2];

# %%
col1 = foobar[Int(c1)]

# %%
@btime col1[2];

# %%
foobar

# %%
foobar[2]

# %%
eltype(foobar)

# %%
tt = Vector{Vector{T}} where T

# %%
function row(arr::Vector{Vector{T}} where {T}, i::Int64)
    res = Array{Any}(undef, length(arr))
    for c in eachindex(arr)
        @inbounds res[c] = @inbounds arr[c][i]
    end
    res
end

function row2(arr::Vector{Vector{T}} where {T}, i::Int64)
    # [r[i] for r in arr]
    map(r -> @inbounds(r[i]), arr)
end

function row3(arr, i)
    res = eltype(arr)(undef, length(arr))
    n = 1
    for j in arr
        @inbounds res[n] = @inbounds(j[i]) 
        n += 1
    end
    return res
end

function row4(arr::Vector{Vector{T}} where {T}, i::Int64)
    [@inbounds(r[i]) for r in arr]
    # map(r -> @inbounds(getindex(r, i)), arr)
end

function row5(arr::Vector{Vector{T}} where {T}, i::Int64)
    # [r[i] for r in arr]
    map(n -> @inbounds(arr[n][i]), 1:size(arr,1))
end

# %%
row5(foobar, 2)

# %%
eltype(eltype(foobar))

# %%
@btime row(foobar, 2)

# %%
@btime row2(foobar, 2)

# %%
@btime row3(foobar,2);

# %%
@btime row4(foobar,2);

# %%
@btime row5(foobar,2);

# %%
tt = Table(a=foobar[1], b=foobar[2])

# %%
@btime tt.a;

# %%
@btime foobar[Int(c1)];

# %%
@benchmark foobar[Int(c1)][2] = iv setup=(iv=44)

# %%
@benchmark tt.a[2] = iv setup=(iv=44)

# %%
@btime row2(foobar,2);

# %%
@btime tt[2]

# %%
Base.summarysize(tt)

# %%
Base.summarysize(foobar)

# %%
coltta = tt.a;

# %%
colonefoo = foobar[Int(c1)];


# %%
@btime coltta[2];

# %%
@btime colonefoo[2];

# %%
keys(getfield(locdat,:data))

# %%
getfield(locdat, :data)

# %%
