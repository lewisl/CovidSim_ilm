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
    @assert length(keyarr) == length(valuearr) "Length of sourcearr not equal length of targetarr"
    for i in eachindex(keyarr)
        if x == keyarr[i]
            return valuearr[i]
        end
    end
    return nothing
end  # not used as hardwired if-test mapping is WAY faster


function countvec!(resvec::Vector{Int}, sourcevec, intmapper::Function)
    for val in sourcevec
        resvec[intmapper(val)] += 1
    end
end

function countvec!(resvec::Vector{Int}, sourcevec, mapdict::Dict,  intmapper=mapviadict)
    for val in sourcevec
        resvec[intmapper(mapdict, val)] += 1
    end
end


@inline function mapcondition(x::condition) # from enum to ordinal int
    if x == uninfected
        0
    else
        Int(x)-4
    end
end


@inline function mapcondition(x::Int) # from ordinal int to enum
    if 0 <= x <= 4
        if x == 0
            uninfected
        elseif x == 1
            nil
        elseif x == 2
            mild
        elseif x == 3
            sick
        else # last condition 4
            severe
        end
    else
        @assert false "invalid integer for mapping to condition $x"
    end
end

@inline function mapagegrp(x::agegrp) # from enum to int
    Int(x)
end

function mapstatus(x::status)
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


function maptransition(x::Integer) # from integer column to elements of transition vector
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

const vaxdict = Dict(:Pfizer=>1, :Moderna=>2, :JnJ=>3)
const variantdict = Dict(:base => 1, :alpha=>2, :delta=>3, :omicron_ba1=>4, :omicron_ba2=>5, :omicron_ba4_5=>6)


function mapviadict(mapdict, x::Symbol)
    get(mapdict, x) do
        throw(DomainError(x, "Argument must be one of $(keys(mapdict))"))
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


function repeat_join(l1::Vector, l2::Vector)
    len1 = length(l1)
    len2 = length(l2)
    res = Vector{Symbol}(undef, len1*len2)
    for (i1, it1) in enumerate(l1)
        for (i2, it2) in enumerate(l2)
            res[(i1-1)*len2+i2] = Symbol(it1, "_", it2)
        end
    end
    return res
end


function repeat_join(l1::Union{Symbol, String}, l2::Vector)
    repeat_join([l1], l2)
end
