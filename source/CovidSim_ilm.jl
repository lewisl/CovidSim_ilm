# TODO
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # fix all the travel functions to latest APIs
    # put in an inflection measure

    # add transq to ilm and test
    # recovery should not guarantee future immmunity
    # vaccination with 3 vaccines / 1 or 2 shots
    # support variants for transition, spread and vaccination
    # use currying to simplify case APIs
    # extend to two years
    # clean up reports and notebooks to work with latest ilm model: get rid of some...
        # get rid of riskmx
    # should quarantine be special or is it extreme social distancing--with no contacts?
        #= 
        tricky because we only using contacts for outgoing contacts by spreaders.
        we would need to reject contacts by the recipient ALSO--not that hard
        this would help with viral load modeling
        =#
    # implement "rings" to set boundaries on contacts and create spreader events and high-risk communities
        #=
        this changes spread logic a lot
        =#
    # decouple spread process from current specifics for "contact", "touch" and "infect"



__precompile__(true)

module CovidSim_ilm

# required
using DelimitedFiles
using DataStructures
using OrderedCollections
using OrderedCollections: FrozenLittleDict
using DataFrames
using CSV
using Random
using Distributions
using StatsBase
using Printf
using Plots
using PlotThemes
using StatsPlots
using Dates
using YAML
using TypedTables
using OffsetArrays
using SplitApplyCombine
using Interpolations


# order matters for these includes!
include("dec_tree.jl")
include("setup.jl")
include("tracking.jl")
include("cases.jl")
include("test_and_trace.jl")
include("transition.jl")
include("spread.jl")
include("vax.jl")
include("sim.jl")
include("johns_hopkins_data.jl")

# functions for simulation
export    
    setup,              
    run_a_sim,
    day_ctr,
    isolate!,
    unisolate!,
    grab,
    input!,
    plus!,
    minus!,
    r0_sim,
    set_by_level

# functions for spreading
export
    Infectparams,
    spread!,
    seed!,
    numcontacts,
    istouched,
    isinfected

# functions for vaccines
export
    Vaccineparams,
    Vaxsched,
    makevaxfn,
    getashot!

# functions for transition
export
    transition!,
    dotransition!

# functions for cases
export
    test_and_trace,     
    Spreadcase,
    sd_gen,
    seed_case_gen,
    t_n_t_case_gen,
    case_setter,
    bayes,
    tntq,
    shifter

# functions for setup
export                  
    build_data,
    setup

# functions for tracking
export                  
    reviewdays,
    showq,
    cumplot,
    newplot,
    dayplot,
    dayanimate2,
    review_history,
    make_series,
    virus_outcome

# queues and caches (variables) for tracking
export            
    travelq,
    spreadq,
    transq,
    day2df,
    map2series

# functions for decision trees
export                  
    setup_dt,
    display_tree,
    sanitycheck,
    getseqs

# functions for accessing data from Johns Hopkins
export                 
    get_real_data,
    loc2df,
    read_actual

# control constants
export                  
    age_dist,
    sickdays,
    sickdaylim

# constants for geo data
export      
    fips,
    state,
    size_cat,
    popsize,
    major,
    large,
    medium,
    small,
    smaller,
    rural

# constants for indices to population matrix
export     
    status,         
    unexposed,
    infectious,
    recovered,
    dead,
    notsick,
    condition,
    nil,
    mild,
    sick,
    severe,
    totinfected,
    statuses,
    conditions,
    allconds,
    condnames,
    infectious_cases,
    transition_cases,
    map2series,
    series_colnames,
    agegrp,
    age0_19,
    age20_39,
    age40_59, 
    age60_79, 
    age80_up, 
    agegrps,
    n_agegrps,
    recvrisk,
    totalcol,
    symtocond,
    symtostat,
    symtoage,
    symtoallconds


###########################################################################
# module global constants (except in Julia things aren't really constant!)
###########################################################################


"""
- use incr!(day_ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(day_ctr, :day) to remove :day and return its current value, set it to 0
- use day_ctr[:day] to return current value of day
"""
const day_ctr = counter(Symbol) # from package DataStructures


################################################################
# constants for data structure indices
################################################################

# control constants
const age_dist = [0.251, 0.271, 0.255, 0.184, 0.039]
const sickdaylim = 25
const sickdays = 1:sickdaylim   # rows

# geo data: fips,county,city,state,sizecat,pop,density
const fips = 1
const county = 2
const city = 3
const state = 4
const sizecat = 5
const popsize = 6
const density = 7
const anchor = 8
const restrict = 9
const density_fac = 10

# population centers sizecats
const major = 1  # 20
const large = 2  # 50
const medium = 3
const small = 4
const smaller = 5
const rural = 6

# enum values for condition, status and agegrp to use in popmatrix
@enum condition begin
    notsick=0 
    nil=5 
    mild 
    sick 
    severe
end

@enum status begin
    unexposed=1 
    infectious 
    recovered 
    dead
end

@enum agegrp begin
    age0_19=1 
    age20_39 
    age40_59 
    age60_79 
    age80_up
end

const statuses = collect(instances(status))
const infectious_cases = [nil, mild, sick, severe]
const transition_cases = [recovered, nil, mild, sick, severe, dead]
const allconds = vcat(infectious_cases, statuses) # note excludes notsick::condition=0
const agegrps = instances(agegrp) # tuple of enums
const n_agegrps = length(instances(agegrp))

#= 
lookup tables for enum values: 
- don't need lookup for Int or Symbol: just use Symbol(nil) and Int(nil)-->these are faster than any lookup
- for symbol use symcond[:nil] => nil::condition = 5
- for string use symcond[Symbol("nil")] => nil::condition = 5
=#



inst_c = instances(condition)
"""
    symtocond[cond::Symbol]
Dict used as lookup table to convert symbol or string to enum value for a condition.

Examples:
- for symbol use symtocond[:nil] => nil::condition = 5
- for string use symtocond[Symbol("nil")] => nil::condition = 5
    
"""
const symtocond = freeze(Dict(zip(Symbol.(inst_c), inst_c)))

# lookup table for status enum values
inst_s = instances(status)
"""
    symtostat[stat::Symbol]
Dict used as lookup table to convert symbol or string to enum value for a status.

Examples:
- for symbol use symtostat[:infectious] => infectious::status = 2
- for string use symtostat[Symbol("infectious")] => infectious::status = 2
    
"""
const symtostat = freeze(Dict(zip(Symbol.(inst_s), inst_s)))

# lookup table for agegrp
inst_a = instances(agegrp)

"""
    symtoage[age::Symbol]
Dict used as lookup table to convert symbol or string to enum value for an age.

Examples:
- for symbol use symtoage[:age0_19] => agegrp::age0_19 = 1
- for string use symtoage[Symbol("age0_19")] => agegrp::age0_19 = 1
    
"""
const symtoage = freeze(Dict(zip(Symbol.(inst_a), inst_a))) # .5x time of regular dict

# lookup table for combined status and condition
const symtoallconds = merge(symtostat, symtocond)

# other columns used only in series dataframes
const totinfected       = 9
const travelers         = 10
const isolated          = 11

# stats series/dataframe columns
# columns of history series: first 5 cols are agegrps, 6th is total
const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24, 
                    nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54)
const condnames  = Dict(:unexposed=>"unexposed", :infectious=>"infectious", :recovered=>"recovered", :dead=>"dead",
                    :nil=>"nil", :mild=>"mild", :sick=>"sick", :severe=>"severe", 9=>"totinfected")
const totalcol = 6


# traveling constants
const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group

end # module CovidSim
