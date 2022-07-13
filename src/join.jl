function repeat_join(l1, l2)
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