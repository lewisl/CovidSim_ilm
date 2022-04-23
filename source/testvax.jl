##
using CovidSim_ilm
##
using StatsBase
using TypedTables
using BenchmarkTools
using Distributions
using YAML
using PrettyPrint
using Plots
using SplitApplyCombine
##
cd(joinpath(homedir(),"Dropbox/Covid Modeling/Covid-ILM/source"))
##
seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :base, agegrps)
##
ndays = 100
locale = 38015
##
result_dict, series = run_a_sim(ndays, locale; 
    dovax=false, 
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantsfilename = "variants.yml",
    showr0=false, 
    silent=true, 
    runcases=[seed_1_6]);
##
locdat = result_dict[:dat]["popdat"][38015]
##
any(last.(locdat.vaxrcvd) .!= :none)
##
cumplot(series, locale)
##
any(series[38015][:cum][:, map2series[:pfizer][6]] .!= 0)
##