using PrettyPrint

function dict_key_to_symbol(d)
    Dict(Symbol(k)=>
            (!(typeof(v) <: AbstractDict) ? v : dict_key_to_symbol(v))
        for (k,v) in d)
end

function change_key_type(d; f=Symbol)
    Dict(f(k) =>
        if !(typeof(v) <: AbstractDict)
            v 
        else 
            change_key_type(v, f=f)
        end
        for (k, v) in d)
end

function dict_key_to_string_v2(d)
    Dict(string(k) =>
        if !(typeof(v) <: AbstractDict)
            v
        else
            dict_key_to_string_v2(v)
        end
        for (k, v) in d)
end

# a somewhat deeply nested, ragged dict
d3lvl = Dict("l1_a"=>
                    Dict("l2_a"=>
                        Dict("l3_a"=>5)), 
             "l1_b"=>4.0,
             "l1_c"=>
                    Dict("l2_a"=>
                        Dict("l3_a"=>5.0, "l3_b"=>7.0)))

# using method dispatch
dict_key_to_symbol_v3(d::AbstractDict) =
           Dict(Symbol(k) => dict_key_to_symbol_v3(v) for (k,v) in d)

@inline dict_key_to_symbol_v3(v) = v

# a more pure recursive version, not better--just more pure
# actually, a lot worse

function more_pure_to_symbol(orig_d, new_d=Dict{Symbol, Any}())
    # if isnothing(iterate(orig_d))
    if isempty(orig_d)
        new_d
    else
        pair = pop!(orig_d)
        k = first(pair)
        v = last(pair)
        new_d[Symbol(k)] =      # have to update new_d in place, which is not very "functional"
            if !(typeof(v) <: AbstractDict)
                v
            else
                more_pure_to_symbol(v)
            end
        more_pure_to_symbol(orig_d, new_d) # finish looping across same level keys in the original dict
    end
end

function change_key_type_recursive(orig_d; T, f, new_d=Dict{T, Any}()) 
    step = iterate(orig_d)  # step is (next_element, iterator_state)
    if isnothing(step)
        new_d
    else
        k = step[1][1]   # step[1] will be a pair for the first element of dict orig_d
        v = step[1][2]   # state of the iterator
        new_d[f(k)] =      
            if typeof(v) <: AbstractDict
                change_key_type_recursive(v, T=T, f=f)
            else
                v
            end
        change_key_type_recursive(Base.rest(orig_d, step[2]), T=T, f=f, new_d=new_d)
    end
end 

# tricky way suggested in forum
# uses function broadcasting and 2 methods
change_dict_type(notdict, args...) = notdict

function change_dict_type(d::AbstractDict, ::Type{T}) where {T}
    Dict{T,Any}(
        T.(keys(d)) .=>
            change_dict_type.(values(d), Ref(T))
    )
end

# uses dict comprehension and 2 methods
dict_key_to_symbol_v3(d::AbstractDict, f) =
    Dict(f(k) => dict_key_to_symbol_v3(v, f) for (k, v) in d)

@inline dict_key_to_symbol_v3(v, f) = v