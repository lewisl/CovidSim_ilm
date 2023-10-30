"""These methods of define_trait! provide inputs to build the primary
simulation data structure that holds traits and outcomes for all 
individuals in the simulation.
"""

####################################################################
# for population data
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
# for history series
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