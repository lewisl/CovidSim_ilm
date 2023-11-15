#####################################################################################
# Data mapping: generally for enums to ordinal integers
#       or ordinal integers to enums
#
# If/elseif is the fastest way to do this with only a handful of items; 
#    use a struct as the mapper would be a close second: 10% slower
#    Dictionaries.jl would be close third: 15% slower   ("Dictionary" is the constructor)
#    both of the latter are easier to maintain and document the assignment
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
@inline function mapit(x, keyarr, valuearr)
    @assert length(keyarr) == length(valuearr) "Length of sourcearr not equal length of targetarr"
    for i in eachindex(keyarr)
        if x == keyarr[i]
            return valuearr[i]
        end
    end
    return nothing
end  # not used as hardwired if-test mapping is WAY faster


@inline function countvec!(resvec::Vector{Int}, sourcevec, intmapper::Function)
    for val in sourcevec
        resvec[intmapper(val)] += 1
    end
end

@inline function countvec!(resvec::Vector{Int}, sourcevec, mapdict::Dict,  intmapper=mapwithdict)
    for val in sourcevec
        resvec[intmapper(mapdict, val)] += 1
    end
end

#= found in the module definition file CovidSim_ilm.jl
    const TOUCHES = [:unexposed, :recovered, :nil, :mild, :sick, :severe]
    const STATUSES = [:unexposed, :infectious, :recovered, :dead]
    const INFECTIOUS_CASES = [:nil, :mild, :sick, :severe]
    const TRANSITION_CASES = [:recovered, :nil, :mild, :sick, :severe, :dead]
    const AGEGRPS = [:age0_19, :age20_39, :age40_59, :age60_79, :age80_up]
=#

@inline function mapcondition(x::Symbol)  # symbol to int
    findit(x, INFECTIOUS_CASES)   # enclosure is const so performance is good
end

@inline function mapcondition(x::Int)   # int to symbol
    findit(x, INFECTIOUS_CASES)
end

@inline function mapagegrp(x::Symbol)
    findit(x,  AGEGRPS)
end

@inline function mapagegrp(x::Int)
    findit(x, AGEGRPS)
end

@inline function mapstatus(x::Symbol)
    findit(x, STATUSES)
end

@inline function map_progression(x::Symbol)
    findit(x, TRANSITION_CASES)
end

@inline function map_progression(x::Int)
    findit(x, TRANSITION_CASES)
end

@inline function maptouch(x::Symbol)
    findit(x, TOUCHES)
end

@inline function maptouch(x::Int)
    findit(x, TOUCHES)
end

@inline function findit(item::Symbol, vec::Vector{Symbol})  # symbol to int
    for i in eachindex(vec)  # linear search with very few items
        if item == vec[i]
            return i
        end
    end
    return 0
end

@inline function findit(item::Int, vec::Vector{Symbol})   # int to symbol
    return vec[item]
end


const vaxdict = Dict(:Pfizer=>1, :Moderna=>2, :JnJ=>3)
const variantdict = Dict(:base => 1, :alpha=>2, :delta=>3, :omicron_ba1=>4, :omicron_ba2=>5, :omicron_ba4_5=>6)


function mapwithdict(mapdict, x::Symbol)
    get(mapdict, x) do
        throw(DomainError(x, "Argument must be one of $(keys(mapdict))"))
    end
end


function tup2vec(maptup, vals)
    [getfield(maptup, x) for x in vals]
end


"""
    symboltoagegrp(x::Union{Symbol, String})
Lookup a string or symbol that matches an enum value of Enum agegrp.
Generates an error if the string or symbol does not match.

Examples:
- symbol2agegrp(:age0_19) result:  agegrp::age0_19 = 1
- symbol2agegrp("age0_19") result: agegrp::age0_19 = 1
    
"""
function symbol2agegrp(x::Union{Symbol, String})::Symbol    # Agegrp
    x = Symbol(x)  # YAML loads age0_19 as string "age0_19"--> convert to Symbol
    return x
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


# not using this:  close second for performance
# struct conditions
#     nil::Int64
#     mild::Int64
#     sick::Int64
#     severe::Int64
# end
# 
# const mapconds = conditions(1,2,3,4)

# @inline function mapcondition_str(x::Symbol, mapconds=mapconds)::Int64
#     getfield(mapconds,x)
# end