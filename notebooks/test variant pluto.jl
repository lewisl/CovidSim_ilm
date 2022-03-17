### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 9d774e4e-57ea-41ca-96ab-1e67649f3571
push!(LOAD_PATH, joinpath(homedir(), "Dropbox/Covid Modeling/Covid-ILM/source"))

# ╔═╡ 0f10581a-02a6-4fb6-95df-9451149dc39f
using CovidSim_ilm


# ╔═╡ 84663908-688b-4f5c-a260-654b3d94b876
cs = CovidSim_ilm

# ╔═╡ 2b8c21bf-9f04-4066-b0ca-f504b52a3410
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

# ╔═╡ 0b4f1675-ed1b-493e-9a74-342d3ac963ae
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))

# ╔═╡ e9ef00a2-dba9-4f1e-be0c-2a393aedb8d1
md"""
# Test setup and population matrix
"""

# ╔═╡ 1d64b481-c661-4de1-9b54-8ec905c63b37
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

# ╔═╡ fe1afa65-cfe2-4846-9a89-cff0d7bc9228
keys(model)

# ╔═╡ 18618970-1daa-4eaf-8534-9e328b360b36
md"""
### Data tables
"""

# ╔═╡ 3c2e0c5d-d2d4-4e74-87c9-a98cb2d7865a
locdat = model.dat["popdat"][locale]

# ╔═╡ 60734094-307b-4937-9c4f-f4ec5c644248
ages = model.dat["agegrp_idx"][locale]

# ╔═╡ 5971199c-0fe7-48fc-b523-09bad0120f50
columnnames(locdat)

# ╔═╡ be3ad807-725a-476b-bd91-e65838e9f844
countmap(locdat.agegrp)

# ╔═╡ ca7d387b-c048-4d88-916f-bbf490012f07
countmap(locdat.status)  # everyone begins as unexposed

# ╔═╡ 3288a611-7686-45ce-b079-88b3f1eeab6b
md"""
### Social Parameters
"""

# ╔═╡ 8fa4c282-74ae-4f08-bcd5-5af0e4b20a96
fieldnames(typeof(model.social))

# ╔═╡ 3c6b339a-f5e3-42f5-ae32-705442ff99b2
typeof(model.social.gammashape)

# ╔═╡ a1c57cf9-d733-4f60-8298-ed5df7d23d4f
model.social.contactfactors

# ╔═╡ b65f1288-0410-4df0-b486-e7d31632ffb1
model.social.contactfactors[Int(sick)-4, Int(age0_19)]

# ╔═╡ d939f648-6b97-4474-8698-ee82f4e8c449
touchfactors = model.social.touchfactors

# ╔═╡ 3bc34818-0c91-41d9-86b3-e87e77ed868a
md"""
### Parameters for infection based on variant
"""

# ╔═╡ 908e458f-d2dc-441d-9825-23b570f5882b
fieldnames(typeof(model.infectset[:omicron]))

# ╔═╡ f1c62977-4f0d-4390-ab11-eb7bab14213b
model.infectset

# ╔═╡ c8e250ef-ad6d-4aea-abf5-6e2db29f57b5
pprintln(model.infectset)

# ╔═╡ e33e704b-2a81-4ec2-b2f5-b356cdd4c2f5
pprintln(model.infectset[:base])

# ╔═╡ 2d067b66-0f85-45c3-b5e0-1cd68b9c9db5
sendbase = model.infectset[:base].sendrisk
recvbase = model.infectset[:base].recvrisk

# ╔═╡ 3c4aace8-65f4-4ec0-b1d8-d354c55290e9
newomisend = zeros(25)
newomisend[1:11] = 1.5 .* sendbase[1:11]
newomisend

# ╔═╡ aecdf4d1-658a-44cb-a83d-687101a961e8
pprint(model.infectset[:alpha])

# ╔═╡ e4833ff3-2041-4809-9bbe-c5fe9cbe7caf
pprintln(model.infectset[:omicron])

# ╔═╡ 842f51f9-9c0d-4159-aa4b-4a7349466554
# is shifter working?
shifter(touchfactors, (.18, .3)...)[:, Int(age40_59)]

# ╔═╡ bda5e583-d954-4339-aa56-9c5fa18afd58
model.transitionset # the decision transition matrices for all age groups are loaded

# ╔═╡ adb8fbdd-d0f0-4a0c-85f0-6be192ffb7c4
basetransition = model.transitionset[:base]

# ╔═╡ bec9c36f-c6da-4a39-8f91-58568e41a2c0
fieldnames(typeof(basetransition))

# ╔═╡ 5360f08d-9161-4aee-9464-e4dc67c5dfb8
fieldnames(typeof(basetransition.tree))

# ╔═╡ 76aa8f66-905f-4931-913c-2cdbca3943f9
keys(age80tree)  # array of Transitiondef

# ╔═╡ 40670f9f-08c8-4795-a05e-c4a1baca4702
basetransition.factors

# ╔═╡ 60218f76-8a50-4c8c-8b76-5c17074230c7
basetransition.factors.vaxhalflifeadjust

# ╔═╡ bef6aff6-c0e8-40ed-8df2-979ef2ec0181
println(typeof(age80tree))
breakday_idx = 5
println("for breakday $breakday_idx")
println("fieldnames: ", fieldnames(typeof(age80tree[breakday_idx])))
println("sickday ",age80tree[breakday_idx].sickday)
println("transition \n", age80tree[breakday_idx].transition)

# ╔═╡ 99391713-d07d-46e2-9d35-09b4f36de00f
md"""
### Sanity check the transition tree
"""

# ╔═╡ 3668f4ad-6474-4ffb-bd47-c852d65ebdd5
cs.sanitycheck(basetransition.tree)

# ╔═╡ 31d47cb1-1dd4-47ea-8e02-136744b42737
md"""
# Run a simulation
"""

# ╔═╡ c99538d2-300a-49e9-97cb-b5f92275ba76
md"""
### Build the simulation model
"""

# ╔═╡ 355b57aa-294d-4e83-9f18-e97f8baeff19
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

# ╔═╡ 0fc5eaac-4867-4a3a-9286-93dbb6318625
md"""
### Create a seed case
"""

# ╔═╡ cde6aa53-06e1-4e8f-a991-70bebd78938d
seed20_39_day1 = makesickseedfunc(; cond=nil, variant=:base, sickday=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=1, forstartofday=true)
seed40_59_day1 = makesickseedfunc(; cond=nil, variant=:base, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=1, forstartofday=true)                            

# ╔═╡ b3ee39c7-eca3-4022-b9ae-446788011d71
seed20_39_omicron = makesickseedfunc(; cond=nil, variant=:omicron, sickday=1, filter=[Term(:agegrp, age20_39), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=360, forstartofday=true)
seed40_59_omicron = makesickseedfunc(; cond=nil, variant=:omicron, sickday=1, filter=[Term(:agegrp, age40_59), Term(:status, unexposed)], 
                            cnt=3, forlocale=0, forday=360, forstartofday=true)                            

# ╔═╡ 9f76bfd7-7a66-4a93-ad01-f8828add6602
md"""
### Run the simulation model
"""

# ╔═╡ ebde51b1-3793-4d41-bbb4-ded5ae838f64
md"""
### Plot results
"""

# ╔═╡ d75f663f-9298-4f38-8532-3ec8b1acfe13
md"""
Note that the orange line labeled Infectious, which shows the current number of infected people, is *not* what you see in newspaper accounts. In this plot Infectious shows the net infected people: There were some sick people as of the day before. Some more people got sick today. Some people got better: they're not infectious any more--they recovered and are on the blue line. Sadly, some people died--they're not infectious either--they're dead and are on the green line. So net infected is yesterday + new today - recovered today - died today. Newspaper tracking shows the new infections of each day--who got sick today? Tomorrow, if no one new got sick the line would be at zero--even though the people who got sick yesterday aren't better yet. So, the newspaper line goes up and down faster. Yet another approach is to show the cumulative number of infected people: This keeps going up until no one new gets infected--then the line is high but levels off. 
"""

# ╔═╡ 94df173f-afc2-4c43-ab2e-10e22c9d8fbf
(cs.vaxlist, cs.variantlist)

# ╔═╡ 3b93689c-ad50-458d-a299-dfde6bab573e
keys(model)

# ╔═╡ 6fc90702-afab-4659-9686-f96a5b709bb0
model.dat

# ╔═╡ 327b7929-b734-4deb-9bb5-ec66fd6abbea
locdat = popdat[locale]
countmap(locdat.cond)

# ╔═╡ c944828e-2d60-419d-8940-eaf7250c4466
countmap(locdat.status)

# ╔═╡ d40c6633-3de7-4021-898b-35ba33abce89
countmap(locdat.variant)

# ╔═╡ 917568b0-2dbf-4c82-b0ac-add7375fc6be
println(":cum ", typeof(series[locale].cum))
println(":cum ", size(series[locale].cum))
println(":new ", size(series[locale].new))

# ╔═╡ d55902b7-78be-48d8-8ce6-2ef2b7b74034
md"""
## Test a social distancing case
"""

# ╔═╡ 7bf10de8-0d51-4e49-8c24-98f3b03eeaf1
sd1 = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, include_ages=[])    

# ╔═╡ a9bf4ad4-5578-4e9d-8510-5dede3e5da83
sd1_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, include_ages=[])

# ╔═╡ 259c8213-93b2-407b-96d8-d684460f800f
outdat = result_dict.dat["popdat"][locale]
all(outdat.sdcomply .== :none)

# ╔═╡ bd13e708-ead2-42fb-b623-7c81f4381b84
md"""
## Social distancing only among those age40_59, age60_79, age80_plus
"""

# ╔═╡ 90a19908-b839-49fb-ad8f-dd6b11ac9641
sdolder = sd_gen(startday = 55, comply=0.9, cf=(.2,1.0), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# ╔═╡ ad50b993-f8c2-45df-a02c-a4e1c6d908b9
sdolder_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age40_59, age60_79, age80_up])    

# ╔═╡ 225fa53f-3434-497b-a94c-35fa087fa247
olderdat =popdat[locale]

sd = findall(olderdat.sdcomply .!= :none)

@Select(agegrp, sdcomply)(olderdat[sd])

# ╔═╡ 4cf03f24-8506-4d60-9b8e-0fdbcfc4b62d
md"""
## Social Distancing starts with everyone and then the younger folks party
"""

# ╔═╡ bd30645a-88f3-43db-84f6-ad7d9a0d6826
sdyoung_end = sd_gen(startday = 90, comply=0.0, cf=(.2,1.5), tf=(.18,.6), name=:mod_80, 
    include_ages=[age0_19, age20_39])    

# ╔═╡ abc91921-a12f-4808-ae98-04117115bc19
cumplot(series, locale)

# ╔═╡ fb2100ed-a1b4-4de3-8797-6fa8905a6b03
cumplot(series, locale, [:nil, :mild, :sick, :severe])

# ╔═╡ d03454f9-f059-455f-b7b7-2ea1f3cff588
newplot(series, locale)

# ╔═╡ 5a3983ac-d843-4454-91a4-96b552390aba
cumplot(series, locale, [:base, :omicron])

# ╔═╡ 86180e20-8d92-4bc3-8f0f-4ed52e510ad4
colmap = series[locale].cols

# ╔═╡ e94b24a4-5178-48dc-a901-a64d54901079
series[locale].cum[end, colmap[:base]]

# ╔═╡ e05adb58-9196-4ed9-8283-cf392992a46d
series[locale].cum[end, colmap[:base]]

# ╔═╡ ec55190a-a09f-4314-8cb4-f307bf852d1d
keys(popdat)

# ╔═╡ a9b17219-c49e-4d20-9d21-c483dc9381b3
keys(series)

# ╔═╡ 668f274e-3b39-4fb3-9aca-76e3b3cdc0e3
popdat[locale]

# ╔═╡ e7f31288-1be3-45cb-b96e-658c1a53acad
virus_outcome(series, locale, base=:pop)

# ╔═╡ 98f46adc-3057-4285-bcf0-feaeda244df0
fieldnames(typeof(series[locale]))

# ╔═╡ e5bbf6af-2fee-4ade-89f9-47743878327b
map2series = series[locale].cols

# ╔═╡ 10637ffc-6192-4feb-87a8-78cebed66064
series[locale].new[:, map2series[:infectious]]

# ╔═╡ 55989787-8214-4da8-a6a5-675fe59702b1
virus_outcome(series, locale, base=:pop)

# ╔═╡ 83002fad-48ae-4432-9f97-03c72e407882
cumplot(series, locale)

# ╔═╡ e8b9a20c-4b8c-4ee5-834a-a11492068182
cumplot(series, locale,[:infectious, :dead])

# ╔═╡ 9a1cf724-a892-4e16-bfd7-927ef672d903
cumplot(series, locale)

# ╔═╡ 38eca2a8-5b8d-4a67-9eda-b07d6c993c5b
cumplot(series, locale, [:infectious, :dead])

# ╔═╡ e26a0bca-4def-4de9-928d-5029d418a9c1
mixdat = popdat[locale]

sd = findall(mixdat.sdcomply .!= :none)
young = findall((mixdat.agegrp .== age0_19) .| (mixdat.agegrp .== age20_39))
sd_young_idx = intersect(sd, young)
old = findall((mixdat.agegrp .== age40_59) .| (mixdat.agegrp .== age60_79) .| (mixdat.agegrp .== age80_up));

# ╔═╡ bae2ac29-fcb9-4768-8e1a-65963b4ebede
youngtab = Table(mixdat[young])
count(youngtab.sdcomply .== :none)

# ╔═╡ 665400ce-6858-4113-b672-df264f9617b3
oldtab = Table(mixdat[old])
count(oldtab.sdcomply .!= :none)

# ╔═╡ 5344ac99-4fc3-4f41-b296-5e1b8b878384
mixdat = popdat[locale]

# ╔═╡ a7554cd1-a67c-46f5-81bd-f69bd670dfa7
typeof(mixdat.agegrp)

# ╔═╡ 5b6979e4-0625-4670-94ba-a9287936aa8c
cumplot(series, locale, [:infectious, :dead])

# ╔═╡ 2b00d861-b5a9-4bd2-b074-41fa58ed1156
@Select(status, agegrp, cond, sdcomply)(locdat)

# ╔═╡ f0df440e-a2ac-46a7-8ad3-2a7132182b1f
age80tree = basetransition.tree.age80_up[1].sickday

# ╔═╡ 1b375018-f4fc-4945-b818-1a1c28652966
result_dict, series = run_a_sim(ndays, locale, showr0=false, silent=true, runcases=[seed_1_6, sd1, sd1_end]);

# ╔═╡ fb049aea-028a-4d82-987c-a76badd9b7b3
popdat, series = runsim(model;
            dovax=true,
            dovariant = true,
            runcases=[seed20_39_day1, seed40_59_day1, seed20_39_omicron, seed40_59_omicron]   # or seed_1_6 if using old way
            );

# ╔═╡ c131a2e7-efb4-493a-a289-5686761b4c06
age80tree = basetransition.tree.age80_up[1].transition

# ╔═╡ 7b495bc8-0666-4941-80e8-b8bce0bc802d
popdat, series = runsim(model,
    runcases=[seed_1_6, sd1, sdyoung_end]);

# ╔═╡ c9e88fb9-5c0f-4bae-9e64-4d51eade438e
popdat, series = runsim(model; 
                    runcases=[seed_1_6, sdolder, sdolder_end]);

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╠═9d774e4e-57ea-41ca-96ab-1e67649f3571
# ╠═0f10581a-02a6-4fb6-95df-9451149dc39f
# ╠═84663908-688b-4f5c-a260-654b3d94b876
# ╠═2b8c21bf-9f04-4066-b0ca-f504b52a3410
# ╠═0b4f1675-ed1b-493e-9a74-342d3ac963ae
# ╟─e9ef00a2-dba9-4f1e-be0c-2a393aedb8d1
# ╠═1d64b481-c661-4de1-9b54-8ec905c63b37
# ╠═fe1afa65-cfe2-4846-9a89-cff0d7bc9228
# ╟─18618970-1daa-4eaf-8534-9e328b360b36
# ╠═3c2e0c5d-d2d4-4e74-87c9-a98cb2d7865a
# ╠═60734094-307b-4937-9c4f-f4ec5c644248
# ╠═5971199c-0fe7-48fc-b523-09bad0120f50
# ╠═be3ad807-725a-476b-bd91-e65838e9f844
# ╠═ca7d387b-c048-4d88-916f-bbf490012f07
# ╟─3288a611-7686-45ce-b079-88b3f1eeab6b
# ╠═8fa4c282-74ae-4f08-bcd5-5af0e4b20a96
# ╠═3c6b339a-f5e3-42f5-ae32-705442ff99b2
# ╠═a1c57cf9-d733-4f60-8298-ed5df7d23d4f
# ╠═b65f1288-0410-4df0-b486-e7d31632ffb1
# ╠═d939f648-6b97-4474-8698-ee82f4e8c449
# ╟─3bc34818-0c91-41d9-86b3-e87e77ed868a
# ╠═908e458f-d2dc-441d-9825-23b570f5882b
# ╠═f1c62977-4f0d-4390-ab11-eb7bab14213b
# ╠═c8e250ef-ad6d-4aea-abf5-6e2db29f57b5
# ╠═e33e704b-2a81-4ec2-b2f5-b356cdd4c2f5
# ╠═2d067b66-0f85-45c3-b5e0-1cd68b9c9db5
# ╠═3c4aace8-65f4-4ec0-b1d8-d354c55290e9
# ╠═aecdf4d1-658a-44cb-a83d-687101a961e8
# ╠═e4833ff3-2041-4809-9bbe-c5fe9cbe7caf
# ╠═842f51f9-9c0d-4159-aa4b-4a7349466554
# ╠═bda5e583-d954-4339-aa56-9c5fa18afd58
# ╠═adb8fbdd-d0f0-4a0c-85f0-6be192ffb7c4
# ╠═bec9c36f-c6da-4a39-8f91-58568e41a2c0
# ╠═5360f08d-9161-4aee-9464-e4dc67c5dfb8
# ╠═f0df440e-a2ac-46a7-8ad3-2a7132182b1f
# ╠═c131a2e7-efb4-493a-a289-5686761b4c06
# ╠═76aa8f66-905f-4931-913c-2cdbca3943f9
# ╠═40670f9f-08c8-4795-a05e-c4a1baca4702
# ╠═60218f76-8a50-4c8c-8b76-5c17074230c7
# ╠═bef6aff6-c0e8-40ed-8df2-979ef2ec0181
# ╟─99391713-d07d-46e2-9d35-09b4f36de00f
# ╠═3668f4ad-6474-4ffb-bd47-c852d65ebdd5
# ╟─31d47cb1-1dd4-47ea-8e02-136744b42737
# ╟─c99538d2-300a-49e9-97cb-b5f92275ba76
# ╠═355b57aa-294d-4e83-9f18-e97f8baeff19
# ╟─0fc5eaac-4867-4a3a-9286-93dbb6318625
# ╠═cde6aa53-06e1-4e8f-a991-70bebd78938d
# ╠═b3ee39c7-eca3-4022-b9ae-446788011d71
# ╟─9f76bfd7-7a66-4a93-ad01-f8828add6602
# ╠═fb049aea-028a-4d82-987c-a76badd9b7b3
# ╟─ebde51b1-3793-4d41-bbb4-ded5ae838f64
# ╠═abc91921-a12f-4808-ae98-04117115bc19
# ╟─d75f663f-9298-4f38-8532-3ec8b1acfe13
# ╠═fb2100ed-a1b4-4de3-8797-6fa8905a6b03
# ╠═d03454f9-f059-455f-b7b7-2ea1f3cff588
# ╠═5a3983ac-d843-4454-91a4-96b552390aba
# ╠═86180e20-8d92-4bc3-8f0f-4ed52e510ad4
# ╠═e94b24a4-5178-48dc-a901-a64d54901079
# ╠═e05adb58-9196-4ed9-8283-cf392992a46d
# ╠═94df173f-afc2-4c43-ab2e-10e22c9d8fbf
# ╠═3b93689c-ad50-458d-a299-dfde6bab573e
# ╠═ec55190a-a09f-4314-8cb4-f307bf852d1d
# ╠═a9b17219-c49e-4d20-9d21-c483dc9381b3
# ╠═6fc90702-afab-4659-9686-f96a5b709bb0
# ╠═668f274e-3b39-4fb3-9aca-76e3b3cdc0e3
# ╠═327b7929-b734-4deb-9bb5-ec66fd6abbea
# ╠═c944828e-2d60-419d-8940-eaf7250c4466
# ╠═d40c6633-3de7-4021-898b-35ba33abce89
# ╠═e7f31288-1be3-45cb-b96e-658c1a53acad
# ╠═98f46adc-3057-4285-bcf0-feaeda244df0
# ╠═e5bbf6af-2fee-4ade-89f9-47743878327b
# ╠═917568b0-2dbf-4c82-b0ac-add7375fc6be
# ╠═10637ffc-6192-4feb-87a8-78cebed66064
# ╟─d55902b7-78be-48d8-8ce6-2ef2b7b74034
# ╠═7bf10de8-0d51-4e49-8c24-98f3b03eeaf1
# ╠═a9bf4ad4-5578-4e9d-8510-5dede3e5da83
# ╠═1b375018-f4fc-4945-b818-1a1c28652966
# ╠═55989787-8214-4da8-a6a5-675fe59702b1
# ╠═83002fad-48ae-4432-9f97-03c72e407882
# ╠═e8b9a20c-4b8c-4ee5-834a-a11492068182
# ╠═259c8213-93b2-407b-96d8-d684460f800f
# ╟─bd13e708-ead2-42fb-b623-7c81f4381b84
# ╠═90a19908-b839-49fb-ad8f-dd6b11ac9641
# ╠═ad50b993-f8c2-45df-a02c-a4e1c6d908b9
# ╠═c9e88fb9-5c0f-4bae-9e64-4d51eade438e
# ╠═225fa53f-3434-497b-a94c-35fa087fa247
# ╠═9a1cf724-a892-4e16-bfd7-927ef672d903
# ╠═38eca2a8-5b8d-4a67-9eda-b07d6c993c5b
# ╟─4cf03f24-8506-4d60-9b8e-0fdbcfc4b62d
# ╠═bd30645a-88f3-43db-84f6-ad7d9a0d6826
# ╠═7b495bc8-0666-4941-80e8-b8bce0bc802d
# ╠═e26a0bca-4def-4de9-928d-5029d418a9c1
# ╠═bae2ac29-fcb9-4768-8e1a-65963b4ebede
# ╠═665400ce-6858-4113-b672-df264f9617b3
# ╠═a7554cd1-a67c-46f5-81bd-f69bd670dfa7
# ╠═5344ac99-4fc3-4f41-b296-5e1b8b878384
# ╠═5b6979e4-0625-4670-94ba-a9287936aa8c
# ╠═2b00d861-b5a9-4bd2-b074-41fa58ed1156
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
