using TypedTables

function create_samples(n=10_000)
    println("Creating typedtable")
    status = fill(1,n)
    agegrp = rand(1:5,10000)
    cond = zeros(Int,10000)
    duration = zeros(Int,10000)

    dat = Table(status=status,agegrp=agegrp,cond=cond, duration=duration)
end


function create_changes_tuple(n=5)
    println("Creating vector of changes as tuple")
    changes = Vector{Tuple{Int64, Symbol, Int64}}(undef, 5)
    for i in 1:n
        sym = isodd(i) ? :cond : :duration
        changes[i] = (i*3, sym, i+1 )
    end
    changes
end


function create_changes_named_tuple(n=5)
    println("Creating vector of changes as namedtuple")
    changes = Vector{NamedTuple{(:row, :col, :new), Tuple{Int64, Symbol, Int64}}}(undef,n)
    for i in 1:n
        sym = isodd(i) ? :cond : :duration
        changes[i] = (row = i*3, col=sym, new=i+1 )
    end   
    changes 
end


function map_tuple!(dat, changes)
    map(changes) do change
        getproperty(dat, change[2])[change[1]] = change[3]
    end    
end


function map_named_tuple!(dat, changes)
    map(changes) do change
        getproperty(dat, change.col)[change.row] = change.new
    end    
end


function loop_tuple!(dat, changes)
    for i in 1:length(changes)
        change = changes[i]
        getproperty(dat, change[2])[change[1]] = change[3]
    end         
end


function loop_named_tuple!(dat, changes)
    for i in 1:length(changes)
        change = changes[i]
        getproperty(dat, change.col)[change.row] = change.new
    end          
end


function run_benchmarks()
end