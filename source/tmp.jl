using TypedTables
using PrettyPrint
using DataStructures

const day_ctr = counter(Symbol) # from package DataStructures

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

function bycolumn(num, col1, col2, col3)
    (col1[num], col2[num], col3[num])
end

function bytable(num, tab, sym1, sym2, sym3)
    (getproperty(tab, sym1)[num], getproperty(tab, sym2)[num], getproperty(tab, sym3)[num])
end

function bytup(num, tup, sym1, sym2, sym3)
    (getproperty(tup, sym1)[num], getproperty(tup, sym2)[num], getproperty(tup, sym3)[num])
end

function setcolgroups(tab, colvec::Vector{Symbol})
    colgrp = [getproperty(tab, c) for c in colvec]
    NamedTuple{Tuple(colvec)}(colgrp)
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


function multifilt(dat::T, filt) where T <: Table

    reduce(.&, [getproperty(dat, Symbol(t.trait)) .== t.val for t in filt])
        
end

function setval!(dat::T, filt, totrait) where T <: Table
    # for rowidx in filt
    #     for t in totrait
    #         getproperty(dat, t.trait)[rowidx] = t.val
    #     end
    # end

    for t in totrait
        getproperty(dat, t.trait)[filt] .= t.val
    end
   
end

#=
locdat[(getproperty(locdat, Symbol(status)) .== unexposed) .& (getproperty(locdat, Symbol(agegrp)) .== age20_39)]
=#




            Base.@kwdef struct Term
                trait::Symbol             # column of the population matrix
                val::Union{Enum, Int, Symbol}    # its value for a row
            end

# TODO
       # limit changeable and filterable columns
       # create case statement for follow-on changes
Base.@kwdef struct Seedset
        filter::Vector{Term}
        change::Vector{Term}
        cnt::Int

        function Seedset(filter::Vector{Term}, totrait::Vector{Term}, cnt::Int)
            allowed_filter_columns = [:cond, :status, :agegrp, :variant, :vaxstatus, :tested, :quar]
            allowed_change_columns = [:cond, :variant, :sickday, :recovered, :vaxrcvd, :quar, :tested]
                #=
                Notes: 
                    :cond used to seed simulation with virus in people
                    :cond goes with :variant and :sickday--if not provided,
                        :variant set to :base and :sickday to 1
                    :recovered may be used to test operation of simulation with assumed context

                    NOT IMPLEMENTED YET:
                    :quar and :tested may be used to implement cases
                    :sdcomply may used to implement cases
                =#
            for f in filter
                if in(f.trait, allowed_filter_columns)
                else
                    throw(DomainError(f.trait, "Column not allowed to be filtered.")) 
                end
            end

            for tt in change
                if in(tt.trait, allowed_change_columns)
                else
                    throw(DomainError(tt.trait, "Column not allowed to be changed.")) 
                end
            end

            new(filter, change, cnt)
        end

end


function settrait!(locdat, day, s::Seedset; startofday=startofday)

    startofday || return

    if day == day_ctr[:day]

        # indices that meet filter criteria
        filt = findall(reduce(.&, [getproperty(locdat, Symbol(f.trait)) .== f.val for f in s.filter]))

        if isempty(filt) 
            @warn "No people met filter criteria. No seeding occurred. Proceeding."
            return
        end

        cnt = min(length(filt), s.cnt) # only do as many as requested or as we have

        filt = cnt < length(filt) ? filt[1:cnt] : filt

        for ch in s.change
            getproperty(locdat, ch.trait)[filt] .= ch.val  # set the value for each trait column
        end

    end
end


function dotest(;day = 1, n=100)

    reset!(day_ctr, :day)  # return and reset key to 0 :day leftover from prior runs
    inc!(day_ctr, :day, day)

    # create a population table
    pop = Dict(48000=>pop_data(100))  # like sim--dict of locale::Int => Table

    # create a seed
    s1 = Seedset(filter=[term(:agegrp, age20_39), term(:status, unexposed)],
              totrait=[term(:cond, nil), term(:variant, :base), term(:sickday, 1)],
              cnt = 5
            )

    # see if the filter works
    filt_time = @elapsed filt = multifilt(pop, s1.filter)

    # see if updating values works
    update_time = @elapsed setval!(pop, filt, s1.totrait)

    @show filt_time, update_time
    return pop

end


"""
Generate seeding cases.

1. Create a Seedset: a filter for who is being changed; the traits that are being changed; the number of people to be seeded.
- use Seedset(filter::Vector{Term}, totrait::Vector{Term}, cnt::Int)
- a Term is Term(trait::Symbol, val::Union{Enum, Symbol, Int})
- even if you have only one filter or totrait, you must wrap it in a vector
- Ex:
```

```

2. Generate the case callback function with seed_case_gen:

3. run the simulation by including the case callback function in the runcases parameter:

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(day, ss::Seedset) # these args go into the called/returned seed! case
    # Arguments to seed_case_gen can be arguments to the inner function that will be run for the case.
    # These arguments are likely specific to what the case will actually do during the simulaton.

    # runcase is called with these params in the case running loop. These arguments are basic parameters
    # of the simulation that have already been set. These won't change for the type of case because they are
    # dependent on the design of the simulation.

    function caserunner(locale, dat, socialparams, infectset, sdcases, ages; startofday)  # args must match runcases loop in run_a_sim
        settrait!(dat, locale, day, s::Seedset; startofday=startofday)  # payload: this is what the function will do when run
    end
end






############################
#   old way--for reference
############################

"""
    seed!(day, cnt, sickday, conds, agegrps, locale, dat)

This is the action function that setups and implements a seeding case all in one execution.
"""
function seed!(day, cnt, sickday, conds, variants, agegrps, locdat; startofday)

    startofday || return  # if true->don't return.  if false->then do return

    @assert length(sickday) == 1 "input only one sickday value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == day_ctr[:day]
        println("*** seed day $(day_ctr[:day]): $(sum(cnt)) $conds to $locale")
        # @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
        make_sick!(locdat; cnt=cnt, ages=agegrps, tocond=conds, tovariant=variants, tosickday=sickday)
    end
end


# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, :base, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, :base, agegrps)
