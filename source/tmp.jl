using TypedTables

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

const agegrps = instances(agegrp)

Base.@kwdef struct term
    trait::Symbol
    val::Union{Enum, Int, Symbol}
end

Base.@kwdef struct seed
        filter::Vector{term}
        totrait::Vector{term}
        cnt::Int
end
        
function pop_data(pop; age_dist=age_dist)

    parts = apportion(pop, age_dist)

    # must use comprehension to initialize vector of vector NOT fill--fill creates identical vectors
    dat = Table(
        pid = collect(1:pop),  # ordinal persistent id for persons in matrix
        status = fill(unexposed, pop),    
        agegrp = reduce(vcat,[fill(age, parts[Int(age)]) for age in agegrps]), 
        cond = fill(uninfected, pop),
        sickday = zeros(Int, pop),   
        variant = [[:none] for _ in 1:pop],
        recovday = [[0] for _ in 1:pop],  
        deadday = zeros(Int, pop),   
        cluster = zeros(Int, pop), 
        sdcomply = fill(:none, pop),  
        vaxstatus = fill(:none, pop),  # :none, :first, :multiple, :full, :booster  maybe others later...
        vaxrcvd = [[:none] for _ in 1:pop],    # vaccine symbols  :pfizer, :moderna, :jnj
        vaxday = [[0] for _ in 1:pop], 
        fullvaxday = zeros(Int, pop),
        tested = falses(pop),  
        testday = zeros(Int, pop),  
        quar = falses(pop),
        quarday = zeros(Int, pop))

    return dat       
end 

const age_dist = [0.251, 0.271, 0.255, 0.184, 0.039]

function apportion(x::Int, splits::Array)
    @assert isapprox(sum(splits), 1.0)
    maxidx = argmax(splits)
    parts = round.(Int, splits .* x)
    diff = sum(parts) - x
    parts[maxidx] -= diff
    return parts
end


function multifilt(s::seed, dat::T) where T <: Table

    dat[ 
        reduce(&, [getproperty(dat, t.trait) .== t.val for t in s.filter])
        , :]

end