using TypedTables
using PrettyPrint
using DataStructures
using Random
using StructArrays


#########################################
#  subset of simulation code to enable testing
#########################################

const DAY_CTR = counter(Symbol) # from package DataStructures

const AGE_DIST = [0.251, 0.271, 0.255, 0.184, 0.039]

@enum Condition begin
    uninfected=0 
    nil=5 
    mild   # 6
    sick   # 7
    severe # 8
end

@enum Status begin
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

const AGEGRPS = instances(agegrp)

        
"""
Pre-allocate and initialize population data for one locale in the simulation.
Returns a TypedTable which is a tuple of arrays:
- each column is a trait of people
- rows are days of the simulatoin
"""
function pop_data1(pop; age_dist=AGE_DIST)

        parts = apportion(pop, age_dist)

        # must use comprehension to initialize vector of vector NOT fill--fill creates identical vectors
        dat = Table(
            status = fill(unexposed, pop),                                          # enum status
            agegrp = reduce(vcat,[fill(age, parts[Int(age)]) for age in AGEGRPS]),  # enum agegrp
            cond = fill(uninfected, pop),                                           # enum Condition
            duration = zeros(Int, pop),                                             # Int
            variant = [Symbol[] for _ in 1:pop],                                     # Vector{Symbol}
            sickday = [[0] for _ in 1:pop],                                         # Vector{Int}
            recovday = [[0] for _ in 1:pop],                                        # Vector{Int}
            deadday = zeros(Int, pop),                                              # Int
            ring = zeros(Int, pop),                                                 # Int (not used as yet)
            sdcase = fill(:none, pop),                                            # Symbol
            vaxstatus = fill(:none, pop),          # :none, :first, :full, :booster  maybe others later...
            vaxrcvd = [[:none] for _ in 1:pop],    # Vector{Symbol} of vaccine symbols  :Pfizer, :Moderna, :JnJ
            vaxday = [[0] for _ in 1:pop],                                          # Vector{Int}
            tested = falses(pop),                                                   # Bool
            testday = zeros(Int, pop),                                              # Vector{Int}
            quar = falses(pop),                                                     # Bool
            quarday = zeros(Int, pop))                                              # Int

    return dat       
end


function pop_data2(pop; age_dist=AGE_DIST)

    parts = apportion(pop, age_dist)

    # must use comprehension to initialize vector of vector NOT fill--fill creates identical vectors
    dat = StructArray(
        status = fill(unexposed, pop),                                          # enum status
        agegrp = reduce(vcat,[fill(age, parts[Int(age)]) for age in AGEGRPS]),  # enum agegrp
        cond = fill(uninfected, pop),                                           # enum Condition
        duration = zeros(Int, pop),                                             # Int
        variant = [Symbol[] for _ in 1:pop],                                     # Vector{Symbol}
        sickday = [[0] for _ in 1:pop],                                         # Vector{Int}
        recovday = [[0] for _ in 1:pop],                                        # Vector{Int}
        deadday = zeros(Int, pop),                                              # Int
        ring = zeros(Int, pop),                                                 # Int (not used as yet)
        sdcase = fill(:none, pop),                                            # Symbol
        vaxstatus = fill(:none, pop),          # :none, :first, :full, :booster  maybe others later...
        vaxrcvd = [[:none] for _ in 1:pop],    # Vector{Symbol} of vaccine symbols  :Pfizer, :Moderna, :JnJ
        vaxday = [[0] for _ in 1:pop],                                          # Vector{Int}
        tested = falses(pop),                                                   # Bool
        testday = zeros(Int, pop),                                              # Vector{Int}
        quar = falses(pop),                                                     # Bool
        quarday = zeros(Int, pop))                                              # Int

return dat       
end

function apportion(x::Int, splits::Array)
    @assert isapprox(sum(splits), 1.0)
    maxidx = argmax(splits)
    parts = round.(Int, splits .* x)
    diff = sum(parts) - x
    parts[maxidx] -= diff
    return parts
end


###############################################################
# code finagling
###############################################################


#########################################
# approach to seeding
#########################################


function multifilt(dat::T, filt) where T <: Table

    reduce(.&, [getproperty(dat, Symbol(t.trait)) .== t.val for t in filt])
        
end

function setval!(dat::T, filt, totrait) where T <: Table

    for t in totrait
        getproperty(dat, t.trait)[filt] .= t.val
    end
   
end


#=
Notes: 
    :cond used to seed simulation with virus in people
    :cond goes with :variant and :duration--if not provided,
        :variant set to :base and :duration to 1
    :recovered may be used to test operation of simulation with assumed context

    NOT IMPLEMENTED YET:
    :quar and :tested may be used to implement cases
    :sdcase may used to implement cases


locdat[(getproperty(locdat, Symbol(status)) .== unexposed) .& (getproperty(locdat, Symbol(agegrp)) .== age20_39)]
=#


Base.@kwdef struct Term
    trait::Symbol             # column of the population matrix
    val::Union{Enum, Int, Symbol}    # its value for a row
end

Base.@kwdef struct Cval
    col::Vector
    val::Union{Symbol, Enum, Int, Bool}
end


Base.@kwdef struct Seedset
        # type::Symbol
        filter::Vector{Term}
        change::Vector{Term}
        cnt::Int

        function Seedset(filter::Vector{Term}, change::Vector{Term}, cnt::Int) 
            allowed_filter_columns = [:cond, :status, :agegrp, :variant, :vaxstatus, :vaxrcvd, :tested, :quar]
            allowed_change_columns = [:cond, :status, :duration, :variant, :vaxstatus, :vaxrcvd]

            for f in filter
                if in(f.trait, allowed_filter_columns)
                else
                    throw(DomainError(f.trait, "Column not allowed to be filtered.")) 
                end
            end

            for ch in change
                if in(ch.trait, allowed_change_columns)
                else
                    throw(DomainError(ch.trait, "Column not allowed to be changed.")) 
                end
            end

            new(filter, change, cnt)
        end

end


function settraits!(locdat, s::Seedset)

    # indices that meet filter criteria
    filt = findall(reduce(.&, [getproperty(locdat, f.trait) .== f.val for f in s.filter]))

    if isempty(filt) 
        @warn "No people met filter criteria for seeding. Proceeding."
        return
    end

    filt = s.cnt < length(filt) ? filt[1:s.cnt] : filt

    for ch in s.change
        col = getproperty(locdat, ch.trait)
        if isa(first(col), Vector) # if element of the column is a vector, push! the value
            for idx in filt
                push!(col[idx], ch.val)
            end
        else 
            col[filt] .= ch.val  # set the value for each trait column
        end
    end

end


function dotest(;day = 1, n=100)

    reset!(DAY_CTR, :day)  # return and reset key to 0 :day leftover from prior runs
    inc!(DAY_CTR, :day, day)

    # create a population table
    pop = Dict(48000=>pop_data(100))  # like sim--dict of locale::Int => Table

    # create a seed
    s1 = Seedset(filter=[term(:agegrp, age20_39), term(:status, unexposed)],
              totrait=[term(:cond, nil), term(:variant, :base), term(:duration, 1)],
              cnt = 5
            )

    # see if the filter works
    filt_time = @elapsed filt = multifilt(pop, s1.filter)

    # see if updating values works
    update_time = @elapsed settraits!(pop, filt, s1.totrait)

    @show filt_time, update_time
    return pop

end


"""
    makesickseedset( ; cond=nil, variant=:base, duration=1, filter::Vector{Term}, cnt)

Return a Seedset that contains the filter for whom to make sick, the traits to be set, and the cnt of people to be changed.

Use this as the first input to seed\\_case\\_gen to create a callback function that will be run
in the simulation loop.
"""
function makesickseedset( ; cond=nil, variant=:base, duration=1, filter::Vector{Term}, cnt)
    Seedset(filter=filter, cnt=cnt, 
            change=[Term(:status, infectious), Term(:cond, cond), Term(:variant, variant), Term(:duration, duration)])
end


"""
    makesickseedfunc( ; cond=nil, variant=:base, duration=1, filter::Vector{Term}, cnt, forlocale=0, forday, startofday)

Create a Seedset that contains the filter for whom to make sick, the traits to be set, and the cnt of people to be changed. 

**And** call seed\\_case\\_gen for you to return the callback function that encloses this Seedset.
"""
function makesickseedfunc( ; cond=nil, variant=:base, duration=1, filter::Vector{Term}, cnt, forlocale=0, forday, forstartofday)
    ss = Seedset(filter=filter, cnt=cnt, 
            change=[Term(:status, infectious), Term(:cond, cond), Term(:variant, variant), Term(:duration, duration)])

    seed_case_gen(ss; forlocale=forlocale, forday=forday, startofday=forstartofday)
end


"""
    makenotsickseedset( ; status, filter::Vector{Term}, cnt)

Return a Seedset that contains the filter for whom to change from sick to either recovered or dead, and the cnt of people to be changed.

Use this as the first input to seed\\_case\\_gen to create a callback function that will be run
in the simulation loop.
"""
function makenotsickseedset( ; status, filter::Vector{Term}, cnt)
    if !in(status, [recovered, dead])
        throw(DomainError(status, "Status must be either recovered or dead"))
    end

    Seedset(filter=filter, cnt=cnt,
            change=[Term(:status, status), Term(:cond, uninfected)])
end


"""
    makenotsickseedfunc( ; status, filter::Vector{Term}, cnt, forlocale=0, forday, startofday)

Create a Seedset that contains the filter for whom to change from sick to either recovered or dead, and the cnt of people to be changed. 

**And** call seed\\_case\\_gen for you to return the callback function that encloses this Seedset.
"""
function makenotsickseedfunc( ; status, filter::Vector{Term}, cnt, forlocale=0, forday, forstartofday)
    if !in(status, [recovered, dead])
        throw(DomainError(status, "Status must be either recovered or dead"))
    end
    
    ss = Seedset(filter=filter, cnt=cnt,
            change=[Term(:status, status), Term(:cond, uninfected)])

    seed_case_gen(ss; forlocale=forlocale, forday=forday, startofday=forstartofday)
end

"""
Generate seeding cases.

1. Create a Seedset: a filter for who is being changed; the traits that are being changed; the number of people to be seeded.
- use Seedset(filter::Vector{Term}, totrait::Vector{Term}, cnt::Int)
- a Term is Term(trait::Symbol, val::Union{Enum, Symbol, Int})
- even if you have only one filter or totrait, you must wrap it in a vector
Example:
```
    ss = Seedset(filter=[Term(trait=:status, val=unexposed)],
             change=[Term(trait=:cond, val=nil)],
             cnt=5)

```julia

2. Generate the case callback function with function `seed\\_case\\_gen`.
Example:
```
    c1 = seed_case_gen(1, ss)  # day 1 and the Seedset above.
```julia

3. run the simulation by including the case callback function in the runcases parameter.
Example:
```
    runsim(popdat, series = runsim(model;
                                dovax=false,
                                runcases=[c1]
                                );
```julia


Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(ss::Seedset; forlocale=0, forday, forstartofday) # these args go into the called/returned seed! case

    function caserunner(locdat, socialparams, infectset, sdcases, age_idx_loc; day, startofday, locale)  # args must match runcases loop in run_a_sim
        if (day == forday) & (startofday == forstartofday) 
            if (forlocale == 0) | (forlocale == locale)
                settraits!(locdat, ss::Seedset)
            end
        end
    end

end

#=
case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=DAY_CTR[:day], startofday=false, locale=loc) 
=#


#########################################
# approach to setting values in population table
#########################################



# array of named tuple
function ant()
end

# named tuple
function nt()
end

# array of tuple
function at()
end

#array of struct
function apstruct(r, pairs::Vector{Cval})  
    for p in pairs
        if isa(p.col[r], Vector)
            push!(p.col[r], p.val)
        else
            p.col[r] = p.val
        end
    end
end

function aptup(r, pairs::Vector{Tuple{Vector{T} where T, Any}})
    for p in pairs
        if isa(p[1][r], Vector)
            push!(p[1][r], p[2])
        else
            p[1][r] = p[2]
        end
    end
end


function aparr(r, pairs::Vector{Vector{Any}}) 
    for p in pairs
        if isa(p[1][r], Vector)
            push!(p[1][r], p[2])
        else
            p[1][r] = p[2]
        end
    end
end


#array of pair with col=>val
function ap(r, pairs::Vector)
    for p in pairs
        if isa(p.first[r], Vector)
            push!(p.first[r], p.second)
        else
            p.first[r] = p.second
        end
    end
end

function apsym(dat, r, pairs)
    for (col, val) in pairs
        vec = getproperty(dat, col)
        if isa(vec[r], Vector)
            push!(vec[r], val)
        else
            vec[r] = val
        end
    end  
end

function aprow(dat, r, changes)
    dat[r] = merge(dat[r],changes)  # need to handle columns that require push!
end

function create_cols(dat)
    c_cond = dat.cond
    c_status = dat.status
    c_duration = dat.duration
    c_variant = dat.variant
    (c_cond=c_cond, c_status=c_status, c_duration=c_duration, c_variant=c_variant)
end

function timesets(; cnt=10000)
    dat = pop_data(cnt)

    c_cond = dat.cond
    c_status = dat.status
    c_duration = dat.duration
    c_variant = dat.variant

    tcolset = 0
    tsymset = 0
    tpairset = 0
    tstructset = 0
    ttupset = 0
    tarrset = 0
    trowset = 0

    for i in 1:Int(cnt/10)
        tcolset += @elapsed begin
                        c_cond[i] = nil
                        c_status[i] = infectious
                        c_duration[i] = 1
                        push!(c_variant[i], :base)
                    end
    end

    st = Int(cnt/10)
    for i in st+1:2*st  
        tpairset += @elapsed begin
            pairs = [c_cond=>nil, c_status=>infectious, c_duration=>1, c_variant=>:base]
            ap(i, pairs)
        end
    end

    for i in 2*st+1:3*st
        tsymset += @elapsed begin
            pairs = [:cond=>nil, :status=>infectious, :duration=>1, :variant=>:base]
            apsym(dat, i, pairs)
        end
    end

    for i in 3*st+1:4*st
        tstructset += @elapsed begin
            pairs = [Cval(c_cond, nil), Cval(c_status, infectious),
                     Cval(c_duration, 1), Cval(c_variant, :base)]
            apstruct(i, pairs)
        end
    end

    for i in 3*st+1:4*st
        ttupset += @elapsed begin
            pairs = [(c_cond, nil), (c_status, infectious),
                     (c_duration, 1), (c_variant, :base)]
            aptup(i, pairs)
        end
    end    

    for i in 4*st+1:5*st
        tarrset += @elapsed begin
            pairs = [[c_cond, nil], [c_status, infectious],
                     [c_duration, 1], [c_variant, :base]]
            aparr(i, pairs)
        end
    end   

    for i in 5*st+1:6*st
        trowset += @elapsed begin
            aprow(dat, i, (status=infectious, cond=nil, duration=1))  # need to build the named tuple programmatically
        end
    end

    println(" tcolset=$tcolset\n tpairset=$tpairset\n tsymset=$tsymset")
    println(" tstructset=$tstructset\n ttupset=$ttupset\n tarrset=$tarrset\n trowset=$trowset")
end



##########
# circular subsets
##############

function circ(n=30)
    infect_idx = collect(1:n)
    contactable_idx = collect(1:3*n)

    shuffle!(contactable_idx)
    taken = pos = 0
    mx = length(contactable_idx)

    for spr in infect_idx
        nc = rand(0:6)  # sloppy use of uniform dist.
        pos = taken + 1
        taken = taken + nc
        if taken <= mx
            sel = pos:taken
        else
            sel = Iterators.flatten((pos:mx, 1:(taken - mx)))
            taken = taken-mx
        end
        println("spreader: ", spr, " nc: ", nc, " sel: ", sel)
        for i in sel
            println("    target: ", contactable_idx[i])
        end
    end
end

##################################################################
#   test performance of closure of arrays
##################################################################

function buildit()
    print("\nEnter n: "); n = parse(Int64, (chomp(readline())))
    return rand(n), collect(1:n)
end

function dostuff!(xfl, xint, arr1, arr2)
    arr1 .*= xfl
    arr2 .+= xint
end

function doitem!(xfl, xint, item, arr1, arr2)
    arr1[item] *= xfl
    arr2[item] += xint
end

arrfl, arrint = buildit()

clos!(xfl, xint) = dostuff!(xfl, xint, arrfl, arrint)
clositem!(xfl, xint, item) = doitem!(xfl, xint, item, arrfl, arrint)

#########################################
# compare dict of vectors to named tuple of vectors
#########################################

function make_dict_vect(num, sz, typ)
    outdict = Dict{Symbol, Vector{typ}}()
    for i in 1:num
        outdict[Symbol("c" * repr(i))] = zeros(typ, sz)
    end
    return outdict
end


function make_nt_vect(num, sz, typ)
    tmp = Dict{String, Vector{typ}}()
    for i in 1:num
        tmp["c" * repr(i)] = zeros(typ, sz)
    end
    return (; (Symbol(k) => v for (k,v) in tmp)...)
end

function dowork!(dat::Dict{Symbol, Vector{Int64}}, cols)
    for i in cols
        dat[Symbol("c" * repr(i))] .+= 3
    end
end

function dowork!(dat::NamedTuple, cols)
    for i in cols
        dat[Symbol("c" * repr(i))] .+= 3
    end
end

function wrap(dat, cols)
    dowork!(dat, cols)
end

