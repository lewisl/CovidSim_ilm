# TODO
    # rewrite R0 sim assuming individual spreading and transition
    # add age_dist as optional parameter in geodata
    # rename transition to be progression
    # add still_infected to summary of statuses
    # should we use statuses of reinfected, breakout? equiv to infectious; need to be filtered
        # whenever we filter for infectious; could help history series
    # do a clean report with vaccines, variants,and social distancing
    # redo vxsched filtervec 
    # per agegrp plots
    # more info
        # get fatality rate by age and co-morbidity CDC, Italian NIH
        # by agegroup, hospitalization %, ICU admission %, fatality %
        # UW virology, expansion of deaths by state on log chart
    # rewrite test and trace to fit new population matrix
    # rewrite quarantine to fit new population matrix--think through social distancing
    # fix all the travel functions to latest APIs
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
    # should sendrisk also depend on condition?  OPTIONALLY, but not for COVID
    



__precompile__(true)

module CovidSim_ilm

# required
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
using Interpolations


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

hash(x::Integer) = uint(x)  # speed up dicts that use integers as keys--especially for transition


################################################################
# constants for data structure indices
################################################################

# control constants
const AGE_DIST = [0.251, 0.271, 0.255, 0.184, 0.039]
const DURATIONLIM = 25
const DURATIONS = 1:DURATIONLIM   # rows

#######################################################################
# enum values for Condition, Status and agegrp to use in population table
#    and related constants
#######################################################################

@enum Condition begin
    uninfected=0 
    nil=5 
    mild   # 6
    sick   # 7
    severe # 8
end

@enum Status begin
    unexposed=1 
    infectious 
    recovered 
    dead
end

@enum Agegrp begin
    age0_19=1 
    age20_39 
    age40_59 
    age60_79 
    age80_up
end


const STATUSES = collect(instances(Status))
const INFECTIOUS_CASES = [nil, mild, sick, severe]
const TRANSITION_CASES = [recovered, nil, mild, sick, severe, dead]
const AGEGRPS = instances(Agegrp) # tuple of enums
const AGEGRPVEC = collect(Symbol.(AGEGRPS)) # vector of symbols
const AGENAMES = vcat(AGEGRPVEC, :total)


# order matters for these includes!
include("data_mapping.jl")
include("dec_tree.jl")
include("setup.jl")
include("tracking.jl")
include("cases.jl")
include("test_and_trace.jl")
include("transition.jl")
include("spread.jl")
include("r0_simulation.jl")
include("vax.jl")
include("sim.jl")
include("johns_hopkins_data.jl")
include("serialize.jl")

# functions for simulation
export    
    buildsim,
    runsim,
    setup,              
    DAY_CTR,
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
    isinfected,
    riskadjust

# functions for vaccines
export
    Vaccineparams,
    Vaxsched,
    Vaxinclude,
    makevaxfn,
    vaccinate!

# functions for transition
export
    transition!,
    dotransition!

# functions for cases
export
    test_and_trace,     
    Spreadcase,
    sd_gen,
    Term,
    Seedset,
    seed_case_gen_old,
    seed_case_gen,
    makesickseedfunc,
    makenotsickseedfunc,
    t_n_t_case_gen,
    case_setter,
    bayes,
    shifter

# functions for setup
export                  
    build_data,
    setup

# functions for tracking
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
    DURATIONS,
    DURATIONLIM

# constants for indices to population matrix
export    
    # enum values for Status and Condition
    Status,         
    unexposed,
    infectious,
    recovered,
    dead,
    uninfected,
    Condition,
    nil,
    mild,
    sick,
    severe,
    STATUSES,
    INFECTIOUS_CASES,
    TRANSITION_CASES,
    # enum and enum values for age groups
    Agegrp,
    age0_19,
    age20_39,
    age40_59, 
    age60_79, 
    age80_up, 
    AGEGRPS


end # module CovidSim
