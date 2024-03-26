# TODO
    # log simulation run messages, or print, or both
    # R0 is broken by shifting spread to one person at a time and update-in-place
    # rewrite test and trace to fit new population matrix
    # test social distancing
    # extend Term to include comparison operations: ==, <=, <, >=, >, in, and not
    # provide a function that filters using Term
    # get rid of old seeding approach
    # implement "rings" to set boundaries on contacts and create spreader events and high-risk communities
        #=
            this changes spread logic a lot
        =#
    # add age_dist as optional parameter in geodata
    # add still_infected to summary of statuses  (?)
    # should we use statuses of reinfected, breakout? equiv to infectious; need to be filtered
        # whenever we filter for infectious; could help history series
        # we could create some track_traits just for tracking purposes:  cleaner
    # do a clean report with vaccines, variants,and social distancing
    # redo vxsched filtervec 
    # per agegrp plots: use term to do flexible filters
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # rewrite quarantine to fit new population matrix--think through social distancing
    # fix all the travel functions to latest APIs
    # should quarantine be special or is it extreme social distancing--with no contacts?
        #= 
            tricky because we only using contacts for outgoing contacts by spreaders.
            we would need to reject contacts by the recipient ALSO--not that hard
            this would help with viral load modeling
        =#
    # should sendrisk also depend on condition?  OPTIONALLY, but not for COVID
        #=
            probably worth building in and "zero-ing" it out for Covid
        =#
    



__precompile__(true)

module CovidSim_ilm

# required
using Tables
using DelimitedFiles
using DataStructures
using OrderedCollections
using OrderedCollections: FrozenLittleDict
using CSV
using Random
using Distributions
using StatsBase
using Printf
using PrettyPrint
using Plots
using PlotThemes
using Dates
using YAML
using TypedTables
using LazyTables
using Interpolations
using NamedTupleTools


######################################################################
# Define module constants and new Base methods
######################################################################

###########################################################################
# module constants (except in Julia things aren't really constant!)
###########################################################################


"""
- use incr!(DAY_CTR, :day) for day of the simulation:  creates and adds 1
- use reset!(DAY_CTR, :day) to remove :day and return its current value, set it to 0
- use DAY_CTR[:day] to return current value of day
"""
const DAY_CTR = counter(Symbol) # from package DataStructures

hash(x::Integer) = uint(x)  # speed up dicts that use integers as keys--especially for progression


################################################################
# constants for data structure indices
################################################################

# control constants
const AGE_DIST = [0.251, 0.271, 0.255, 0.184, 0.039]
const DURATIONLIM = 25   # maximum length of illness in days for anyone
const DURATIONS = 1:DURATIONLIM

#######################################################################
# Symbol values for Condition, Status and agegrp to use in population table
#    and related constants
#
# Symbols are more convenient than enums and faster than strings.
#
# See data_mapping.jl for functions that map symbols to integer values
#######################################################################


#symbols for condition
    :uninfected     # maps to 0
    :nil            # ... 1
    :mild           # ... 2
    :sick           # ... 3
    :severe         # ... 4


# symbols for Status
    :unexposed      # 1
    :infectious     # 2
    :recovered      # 3
    :dead           # 4


# symbols for Agegrps
    :age0_19        # 1
    :age20_39       # 2
    :age40_59       # 3
    :age60_79       # 4
    :age80_up       # 5


# These are VERY IMPORTANT constant vectors used in data_mapping.jl
const TOUCHES = [:unexposed, :recovered, :nil, :mild, :sick, :severe]
const STATUSES = [:unexposed, :infectious, :recovered, :dead]
const INFECTIOUS_CASES = [:nil, :mild, :sick, :severe]
const PROGRESSION_CASES = [:recovered, :nil, :mild, :sick, :severe, :dead]
const AGEGRPS = [:age0_19, :age20_39 , :age40_59, :age60_79, :age80_up]
const AGENAMES = vcat(AGEGRPS, :total)


########################################################################
#  file includes
########################################################################

# order matters for these includes!
include("data_mapping.jl")
include("progression_probs.jl")
include("setup.jl")
include("plotting.jl")
include("simstats.jl")
include("cases.jl")
include("test_and_trace.jl")
include("progression.jl")
include("spread.jl")
include("disease_modeling.jl")
include("r0_simulation.jl")
include("vax.jl")
include("sim.jl")
include("johns_hopkins_data.jl")
include("serialize.jl")


##########################################################################################
# exports
##########################################################################################

# functions for simulation
export    
    buildsim,
    runsim,
    setup_model,
    setup_files,
    setup_yaml,              
    DAY_CTR,
    isolate!,
    unisolate!,
    grab,
    input!,
    plus!,
    minus!,
    r0_sim

# functions for spreading
export
    Infectparams,
    spread!,
    seed!,
    numcontacts,
    istouched,
    isinfected,
    riskadjust

# functions for vaccines
export
    Vaccineparams,
    Vaxsched,
    Vaxinclude,
    makevaxfn,
    vaccinate!

# functions for progression
export
    progression!,
    doprogression!

# functions for cases
export
    test_and_trace,     
    SpreadCase,
    sd_gen,
    Term,
    Seedset,
    seed_case_gen_old,
    seed_case_gen,
    maketraitseedfunc,
    makenotsickseedfunc,
    t_n_t_case_gen,
    case_setter,
    bayes,
    shifter

# functions for setup
export                  
    build_data,
    setup

# functions for plotting
export                  
    reviewdays,
    cumplot,
    newplot,
    dayplot,
    dayanimate2,
    make_series
    
# queues and caches (variables) for tracking
export       
    travelq,
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
    AGE_DIST,
    DURATIONLIM

# constants for indices to population matrix
export    
    # values for Status, Condition and AGEGRP
    STATUSES,
    INFECTIOUS_CASES,
    PROGRESSION_CASES,
    AGEGRPS,
    TOUCHES

# functions for serialization
export
    series_to_csv,
    popdat_to_csv,
    modelinputs_to_yaml,
    model_to_yaml,
    yaml_to_model


end # module CovidSim
