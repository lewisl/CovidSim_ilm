# TODO
    # look for referencing of globals, especially in loops
    # extend Term to include comparison operations: ==, <=, <, >=, >, in, and not
    # provide a function that filters using Term
    # get rid of old seeding approach
    # create generic condition setting function instead of using literals in spread!
    # add age_dist as optional parameter in geodata
    # add still_infected to summary of statuses
    # should we use statuses of reinfected, breakout? equiv to infectious; need to be filtered
        # whenever we filter for infectious; could help history series
    # do a clean report with vaccines, variants,and social distancing
    # redo vxsched filtervec 
    # per agegrp plots: use term to do flexible filters
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

hash(x::Integer) = uint(x)  # speed up dicts that use integers as keys--especially for progression


################################################################
# constants for data structure indices
################################################################

# control constants
const AGE_DIST = [0.251, 0.271, 0.255, 0.184, 0.039]
const DURATIONLIM = 25
const DURATIONS = 1:DURATIONLIM   # rows

#######################################################################
# Symbol values for Condition, Status and agegrp to use in population table
#    and related constants
#
# See data_mapping.jl for functions to map to integer values
#######################################################################


#symbols for condition
    :uninfected     # 0
    :nil            # 1
    :mild           # 2
    :sick           # 3
    :severe         # 4


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


# const STATUSES = collect(instances(Status))
const STATUSES = [:unexposed, :infectious, :recovered, :dead]
const INFECTIOUS_CASES = [:nil, :mild, :sick, :severe]
const TRANSITION_CASES = [:recovered, :nil, :mild, :sick, :severe, :dead]
const AGEGRPS = [:age0_19, :age20_39 , :age40_59, :age60_79, :age80_up]
const AGEGRPVEC = AGEGRPS
const AGENAMES = vcat(AGEGRPVEC, :total)


#################################################################################
#   structs
#################################################################################

# for spread

Base.@kwdef struct Infectparams
    sendrisk::Vector{Float64}
    recvrisk::Vector{Float64}
    recovery_immunity::Dict{Symbol, Float64}
    immunehalflife::Int64
    basemultiplier::Float64
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        function Infectparams(indict::Dict{Symbol, Any})
            Infectparams(
                sendrisk = indict[:sendrisk],
                recvrisk = indict[:recvrisk],
                recovery_immunity = indict[:recovery_immunity],
                immunehalflife = indict[:immunehalflife],
                basemultiplier = indict[:basemultiplier]
                )
        end


Base.@kwdef struct SocialParams
    gammashape::Float64
    indoor_uplift::Float64
    contactfactors::Matrix{Float64}     
    touchfactors::Matrix{Float64}     
end


Base.@kwdef struct SpreadCase       # Base.@kwdef -> use keyword arguments and defaults in constructor
    name::Symbol
    day::Int
    cfdelta::Tuple{Float64,Float64}  
    tfdelta::Tuple{Float64,Float64}  
    comply::Float64             # compliance fraction
    cfcase::Matrix{Float64}
    tfcase::Matrix{Float64}
end

# for progression through disease conditions

Base.@kwdef struct Agetree
    age0_19::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
    age20_39::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
    age40_59::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
    age60_79::Dict{Int, Matrix{Float64}} =  Dict{Int, Matrix{Float64}}()
    age80_up::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()
end

Base.@kwdef struct ProgressionFactors
    riskadjust::Union{Vector{Float64}, Nothing}
    vaxhalflifeadjust::Union{Dict{Symbol, Float64}, Nothing}
    
        # inner method
        function ProgressionFactors(factordict)
            riskadj = get(factordict, :riskadjust, nothing)
            vaxadj = get(factordict, :vaxhalflifeadjust, nothing)
            vaxadj = if !isnothing(vaxadj)
                        Dict(Symbol(k)=>v for (k,v) in vaxadj)
                     end
            new(riskadj, vaxadj)
        end
end

Base.@kwdef struct ProgressionParams
    tree::Union{Agetree, Nothing}
    factors::ProgressionFactors   # use [] for nothing
end


# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
end


########################################################################
#  file includes
########################################################################

# order matters for these includes!
include("data_mapping.jl")
include("dec_tree.jl")
include("setup.jl")
include("plotting.jl")
include("simstats.jl")
include("cases.jl")
include("test_and_trace.jl")
include("progression.jl")
include("spread.jl")
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
    setup,              
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
    DURATIONS,
    DURATIONLIM

# constants for indices to population matrix
export    
    # values for Status, Condition and AGEGRP
    STATUSES,
    INFECTIOUS_CASES,
    TRANSITION_CASES,
    AGEGRPS


end # module CovidSim
