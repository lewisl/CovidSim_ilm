"""These methods of define_trait! provide inputs to build the primary
simulation data structure that holds traits and outcomes for all 
individuals in the simulation.
"""

####################################################################
# for population data:  traits
####################################################################
@kwdef struct Traitdef  # clearer type definition than a tuple
    allowed::Any
    default::Any
    track_history::Bool
end

all_traits = Dict{Symbol, Traitdef}()
# usage: all_traits[:mytrait].values; all_traits[:mytrait].default, all_traits[:mytrait].track_history


# method for a vector of allowed values
# use for Bool with [true, false]
function define_trait!(all_traits, history_traits, trait_symbol, trait_vec, default_val, track_history=true)
    all_traits[trait_symbol] = Traitdef(allowed=trait_vec, default=default_val, track_history=track_history)
end

# method for a bounded range of scalars
function define_trait!(all_traits, history_traits, trait_symbol, min, max, default_val, track_history=true, len=0)
    val = 
        if len == 0
            range(min, max)
        else
            range(min, max, len)
        end
    all_traits[trait_symbol] = Traitdef(allowed=val, default=default_val, track_history=true) 
end

# method for a vector of some type with allowed vector of values
function define_trait!(all_traits, history_traits, trait_symbol, typ::DataType, trait_vec, default_val=nothing, track_history=true)
    defvalue = 
        if isnothing(default_val)
            Vector{typ}()
        else
            Vector{typ}(default_val)
        end
    all_traits[trait_symbol] = Traitdef(allowed=trait_vec, default=defvalue, track_history=track_history)
end

# method for a vector of some type with bounded range of scalar elements
function define_trait!(all_traits, history_traits, trait_symbol, typ::DataType, min, max, len=0, default_val=nothing, track_history=true)
    val = 
        if len == 0
            range(min, max)
        else
            range(min, max, len)
        end
    defval = 
        if isnothing(default_val)
            Vector{typ}()
        else
            Vector{typ}(default_val)
        end

    all_traits[trait_symbol] = Traitdef(allowed=val, default=defval, track_history=track_history)
end


#########################################################################################
# We need to generate these constants that are used throughout the simulation logic,
#   or replace them with better ways to reference these symbols that are part
#   of the model structure.
#
#########################################################################################

    #=
    # These are VERY IMPORTANT constant vectors used in data_mapping.jl
    const TOUCHES = [:unexposed, :recovered, :nil, :mild, :sick, :severe]
    const STATUSES = [:unexposed, :infectious, :recovered, :dead]
    const INFECTIOUS_CASES = [:nil, :mild, :sick, :severe]
    const PROGRESSION_CASES = [:recovered, :nil, :mild, :sick, :severe, :dead]
    const AGEGRPS = [:age0_19, :age20_39 , :age40_59, :age60_79, :age80_up]
    const AGENAMES = vcat(AGEGRPS, :total)

    others includev vectors vaxlist, variantlist, which are built from parameter inputs
    =#


#########################################################################################
# for history series
#########################################################################################

    #=
    Logic already provided to build series columns and series column names.
    - repeat trait by agegrp and total

    Logic to summarize simulation results at end of each time unit (day) 
    - explicitly summarize results by specific traits: statuses, sickness conditions, vaccinations
      and variants

    Some traits are represented as a history of per person outcomes when represented as 
    a vector of vectors, which provides an event history for each person:
    - variant
    - sickday
    - recovday
    - vaxrcvd
    - vaxday
    This is only feasible for traits that only can have a few events per person during the
    entire simulation.

    =#


#########################################################################################
# functions and logic to spread the disease
#########################################################################################
#=
- contacts based on social factors
- touches based on mechanisms that determine and affect transmission of the disease: by fluids, by air as
  droplets or aeresols, by touch and ingestion. Also, characteristics of the pathogen including where it
  can persist outside of human carriers/victims, how long, how pathogen changes within a human that 
  affect its transmission out of carriers
- isinfected based on how a pathogen enters a victim, the victim's biological traits that affect 
  susceptibility to become infected, and the mechanism of infection, latency of pathogen activity 
  in victim that can result in infection
For a simulation, these considerations depend on clinical information and understanding the 
lifecycle of the pathogen but are represented in the simulation by time duration, [SENDRISK] risk of transmitting 
infection, and severity or other conditions that determine how infectious the infected individual 
is to others, and [RECVRISK] the risk of receiving infection based on characteristics of the target individual

=#



#########################################################################################
# functions to progess the disease
#########################################################################################

#=
Covid Model traits
==================
** status
** agegrp
** cond
** duration
** variant
** sickday
** recovday
** deadday
   ring
** sdcase
*^ vaxstatus
*^ vaxrcvd
*^ vaxday
*^ tested
*^ testday
*^ quar
*^ quarday


Legend
=======
** included and required to use with pre-defined values
   included without defined values and not used in simulation logic
*^ included in logic, but use is optional
=#

########## perf test of function passing

function do_math(math_func, x, y)
    math_func(x, y)
end

function adder(x,y)
    +(x,y)
end

function times(x,y)
    *(x,y)
end

# if we demand a required output type of the passed function:
# this trivial example doesn't change performance but can ensure type stability

function do_math_t(math_func, x, y)
    math_func(x,y)::Int64 
end

# how to override an operator to make dict access easier.   <= also makes a good choice
import Base.<=

<=(dct::Dict{Symbol, T}, idx::Symbol) where T = getindex(dct, idx)

# then you can do mydict<=idx