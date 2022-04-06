#####################################################################################
# Data mapping: generally for enums to ordinal integers
#       or ordinal integers to enums
#
#####################################################################################

"""
    mapit(x, keyarr, valuearr)

For an input value x, find the value from the target array where x 
is in the source array. x must match a value in the source array, which 
are semantically keys. The matching value returns an integer index, which selects
the value from the target array. Conceptually, the sourcearr is like keys to a Dict and
the target array are the values of the Dict. For a small number of items, mapit will generally 
be slower than a dict, and may be faster depending on which key is accessed.

Much of the benefit is from not creating a Dict as a mapping. However, if the Dict only
is created once for many accesses then the advantage is small. There is only a benefit
if the lengths of the keys and values is very small.

Ex:
    x = "three"
    sourcearr = ["one", "two", "three", "four"]
    targetarr = [150, 225, 325, 471]
    mapit(x, sourcearr, targetarr) # returns 325
"""
function mapit(x, keyarr, valuearr)
    # @assert length(sourcearr) == length(targetarr) "Length of sourcearr not equal length of targetarr"
    idx = findfirst(isequal(x), keyarr)
    valuearr[idx]
end


function mapcondition(x::condition) # from enum to int
    if x == uninfected
        0
    else
        Int(x)-4
    end
end


function mapcondition(x::Integer) # from int to enum
    if 0 <= x <= 4
        if x == 0
            0
        elseif x == 1
            nil
        elseif x == 2
            mild
        elseif x == 3
            sick
        else
            severe
        end
    else
        @assert false "invalid integer for mapping to condition $x"
    end
end

function mapagegrp(x::agegrp) # from enum to int
    Int(x)
end


function maptouch(x::Union{condition, status}) # from enum to rows of touch parameters
    if x == unexposed
        1
    elseif x == recovered
        2
    elseif x == nil
        3
    elseif x == mild
        4
    elseif x == sick
        5
    elseif x == severe
        6
    else
        @assert false "invalid index to touchfactors $x"
    end
end


function maptransition(x::Union{condition, status}) # from enum to elements of transition vector
    if x == recovered
        1
    elseif x == nil
        2
    elseif x == mild
        3
    elseif x == sick
        4
    elseif x == severe
        5
    elseif x == dead
        6
    else
        @assert false "invalid index for transition vector $x"
    end
end


function maptransition(x::Integer) # from integer to elements of transition vector
    if x == 1 
        recovered
    elseif x == 2
        nil
    elseif x == 3
        mild
    elseif x == 4
        sick
    elseif x == 5
        severe
    elseif x == 6
        dead
    else
        @assert false "invalid index for transition vector $x"
    end
end

function tup2vec(maptup, vals)
    [getfield(maptup, x) for x in vals]
end


#= 
lookup tables for enum values: 
- don't need lookup for Int or Symbol: just use Symbol(nil) and Int(nil)-->these are faster than any lookup
- for symbol use symcond[:nil] => nil::condition = 5
- for string use symcond[Symbol("nil")] => nil::condition = 5
=#

"""
    symboltoagegrp(x::Union{Symbol, String})
Lookup a string or symbol that matches an enum value of Enum agegrp.
Generates an error if the string or symbol does not match.

Examples:
- symbol2agegrp(:age0_19) result:  agegrp::age0_19 = 1
- symbol2agegrp("age0_19") result: agegrp::age0_19 = 1
    
"""
function symbol2agegrp(x::Union{Symbol, String})::agegrp
    x = Symbol(x)
    inst_a = instances(agegrp)
    symtoage = freeze(Dict(zip(Symbol.(inst_a), inst_a))) # .5x time of regular dict
    @assert in(x, keys(symtoage)) "Error: input symbol $x is not an agegrp value."

    symtoage[x]
end


"""
    symbol2condition(x::Union{Symbol String})::condition  
Lookup a string or symbol that matches an enum value of Enum condition.
Generates an error if the string or symbol does not match.

Examples:
- symbol2cond(:nil)  result: nil::condition = 5
- symbol2cond("nil") result: nil::condition = 5
    
"""
function symbol2condition(x::Union{Symbol, String})::condition  
    x = Symbol(x) 
    inst_cond = instances(condition)
    symtocond = freeze(Dict(zip(Symbol.(inst_cond), inst_cond)))

    symtocond[x]
end


"""
    symbol2status(x::Union{Symbol String})::status  
Lookup a string or symbol that matches an enum value of Enum status.
Generates an error if the string or symbol does not match.

Examples:
- symbol2cond(:recovered)  result: recovered::status = 3
- symbol2cond("recovered") result: recovered::status = 3
    
"""
function symbol2status(x::Union{Symbol, String})::status  
    x = Symbol(x) 
    inst_status = instances(status)
    symtostatus = freeze(Dict(zip(Symbol.(inst_status), inst_status)))

    symtostatus[x]
end


"""
    symbol2allconds(x::Union{Symbol String})::Union{condition, status}  
Lookup a string or symbol that matches an enum value of Enum condition or status.
Generates an error if the string or symbol does not match.

Examples:
- symbol2allconds(:nil)  result: nil::condition = 5
- symbol2allconds("nil") result: nil::condition = 5
- symbol2allconds(:dead) result: dead::status = 4
    
"""
function symbol2allconds(x::Union{Symbol, String})::Union{condition, status}  
    x = Symbol(x) 
    inst_status = instances(status)
    inst_cond = instances(condition)

    symtostatus = freeze(Dict(zip(Symbol.(inst_status), inst_status)))
    symtocond = freeze(Dict(zip(Symbol.(inst_cond), inst_cond)))
    symtoallconds = merge(symtostatus, symtocond)

    symtoallconds[x]
end

