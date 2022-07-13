# TODO
    # redo vxsched filtervec 
    # per agegrp plots
    # rewrite R0 sim assuming individual spreading and transition
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


# using Debugger

######################################################################
# Define module constants and new Base methods
######################################################################

###########################################################################
# module constants (except in Julia things aren't really constant!)
###########################################################################


"""
- use incr!(day_ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(day_ctr, :day) to remove :day and return its current value, set it to 0
- use day_ctr[:day] to return current value of day
"""
const day_ctr = counter(Symbol) # from package DataStructures

hash(x::Integer) = uint(x)  # speed up dicts that use integers as keys--especially for transition


################################################################
# constants for data structure indices
################################################################

# control constants
const age_dist = [0.251, 0.271, 0.255, 0.184, 0.039]
const durationlim = 25
const durations = 1:durationlim   # rows

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

#######################################################################
# enum values for condition, status and agegrp to use in population table
#######################################################################

@enum condition begin
    uninfected=0 
    nil=5 
    mild   # 6
    sick   # 7
    severe # 8
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
const allconds = vcat(infectious_cases, statuses) # note excludes uninfected::condition=0
const agegrps = instances(agegrp) # tuple of enums
const n_agegrps = length(instances(agegrp))

# trait columns used in history series
const seriesgroups = [:unexposed, :infectious, :recovered, :dead,         # status
                      :nil, :mild, :sick, :severe, :totinfected,          # infectious cases
                      :Pfizer, :Moderna, :JnJ, :totvaccinated,            # vaccines
                      :base, :alpha, :delta, :omicron_ba1, :omicron_ba2, :omicron_ba4_5]   # variants

# all because building these symbols on the fly is painfully slow because of string catenation!
const seriesbyage =  Dict(:unexposed => (:unexposed_age0_19, :unexposed_age20_39, :unexposed_age40_59, :unexposed_age60_79, :unexposed_age80_up), 
                          :infectious => (:infectious_age0_19, :infectious_age20_39,:infectious_age40_59, :infectious_age60_79, :infectious_age80_up), 
                          :recovered => (:recovered_age0_19, :recovered_age20_39,:recovered_age40_59, :recovered_age60_79, :recovered_age80_up), 
                          :dead => (:dead_age0_19, :dead_age20_39,:dead_age40_59, :dead_age60_79, :dead_age80_up),         
                          :nil => (:nil_age0_19, :nil_age20_39, :nil_age40_59, :nil_age60_79, :nil_age80_up), 
                          :mild => (:mild_age0_19, :mild_age20_39,:mild_age40_59, :mild_age60_79, :mild_age80_up), 
                          :sick => (:sick_age0_19, :sick_age20_39,:sick_age40_59, :sick_age60_79, :sick_age80_up), 
                          :severe => (:severe_age0_19, :severe_age20_39,:severe_age40_59, :severe_age60_79, :severe_age80_up), 
                          :totinfected => (:totinfected_age0_19, :totinfected_age20_39,:totinfected_age40_59, :totinfected_age60_79, :totinfected_age80_up),         
                          :Pfizer => (:Pfizer_age0_19, :Pfizer_age20_39,:Pfizer_age40_59, :Pfizer_age60_79, :Pfizer_age80_up), 
                          :Moderna => (:Moderna_age0_19, :Moderna_age20_39,:Moderna_age40_59, :Moderna_age60_79, :Moderna_age80_up), 
                          :JnJ => (:JnJ_age0_19, :JnJ_age20_39,:JnJ_age40_59, :JnJ_age60_79, :JnJ_age80_up), 
                          :totvaccinated => (:totvaccinated_age0_19, :totvaccinated_age20_39,:totvaccinated_age40_59, :totvaccinated_age60_79, :totvaccinated_age80_up),            
                          :base => (:base_age0_19, :base_age20_39, :base_age40_59, :base_age60_79, :base_age80_up), 
                          :alpha => (:alpha_age0_19, :alpha_age20_39,:alpha_age40_59, :alpha_age60_79, :alpha_age80_up), 
                          :delta => (:delta_age0_19, :delta_age20_39,:delta_age40_59, :delta_age60_79, :delta_age80_up), 
                          :omicron_ba1 => (:omicron_ba1_age0_19, :omicron_ba1_age20_39,:omicron_ba1_age40_59, :omicron_ba1_age60_79, :omicron_ba1_age80_up), 
                          :omicron_ba2 => (:omicron_ba2_age0_19, :omicron_ba2_age20_39, :omicron_ba2_age40_59, :omicron_ba2_age60_79, :omicron_ba2_age80_up),
                          :omicron_ba4_5 => (:omicron_ba4_5_age0_19, :omicron_ba4_5_age20_39, :omicron_ba4_5_age40_59, :omicron_ba4_5_age60_79, :omicron_ba4_5_age80_up)

                        ) 

const serieslevels = Dict(k1 => Dict(k2 => Symbol(k1, "_", k2) 
                                for k2 in [:total, instances(agegrp)...]) 
                          for k1 in seriesgroups)

# other columns used only in series 
const totinfected       = 9
const travelers         = 10
const isolated          = 11

const totalcol = 6


# traveling constants
const travprobs = [1.0, 2.0, 3.0, 3.0, 0.4] # by age group

# order matters for these includes!
include("data_mapping.jl")
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
include("serialize.jl")

# functions for simulation
export    
    buildsim,
    runsim,
    setup,              
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
    seriesgroups,       
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
    age_dist,
    durations,
    durationlim

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
    # enum values for status and condition
    status,         
    unexposed,
    infectious,
    recovered,
    dead,
    uninfected,
    condition,
    nil,
    mild,
    sick,
    severe,

    totinfected,
    statuses,
    conditions,
    allconds,
    infectious_cases,
    transition_cases,
    # enum and enum values for age groups
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
    symbol2conditon,
    symbol2status,
    symbol2agegrp,
    symbol2allconds

# constants for indices to transition arrays
export
    shift,
    recover,
    improve,
    same,
    worse,
    worseplus,
    die


end # module CovidSim
