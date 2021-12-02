function symcol(cols...)
    if !isa(cols, Tuple)
        Symbol(cols)
    else
        first = true
        res = Symbol()
        for col in cols
            if first
                first = false
                res = Symbol(col)
            else
                res = Symbol(res, "__", col)
            end
        end
        return res
    end
end


function symcol2(cols...)
    if !isa(cols, Tuple)
        Symbol(cols)
    else
        first = true
        res = ""
        for col in cols
            if first
                first = false
                res = string(col)
            else
                res *= "__" * string(col)
            end
        end
        return Symbol(res)
    end
end

function symcol4(cols...)
    # @show typeof(cols)
    if length(cols) == 1
        # @show "got here"
        Symbol(cols...)
    else
        first = true
        local res
        for col in cols
            if first
                first = false
                res = Symbol(col)
            else
                res = (res, Symbol(col))
            end
        end
        # @show res
        Symbol(res...)
    end
end


function symcol5(cols...)
    if length(cols) == 1
        Symbol(cols...)
    else
        io = IOBuffer(; maxsize=50)
        first = true
        for col in cols
            if first
                first = false
                # write(io, string(col))
                write(io, Symbol(col))
            else
                write(io, "__")
                # write(io, string(col))
                write(io, Symbol(col))
            end
        end
        seekstart(io)
        res = Symbol(read(io, String))
        close(io)
        return res
    end
end