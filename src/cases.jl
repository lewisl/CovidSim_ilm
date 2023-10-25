####################################################################
# cases.jl
#       pass these cases to run_a_sim in kwarg runcases as a list
#       cases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  
#       then run_a_sim(geofilename, n_days, locales; runcases=cases)
####################################################################




####################################################################
# seeding cases
####################################################################

"""
    seed!(cnt, duration, conds, agegrps, locale, dat)

This is the action function that setups and implements a seeding case all in one execution.
"""
function seed!(cnt, duration, conds, variants, agegrps, locdat)

    @assert length(duration) == 1 "input only one duration value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    println("*** seed day $(DAY_CTR[:day]): $(sum(cnt)) $conds")
    # @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
    make_sick!(locdat; cnt=cnt, ages=agegrps, tocond=conds, tovariant=variants, toduration=duration)

end


# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, :base, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, :base, agegrps)


####################################################################
# new approach to seed people as sick or notsick
#   allows very general filtering
####################################################################

Base.@kwdef struct Term
    trait::Symbol             # column of the population matrix
    val::Union{Enum, Int, Symbol}    # its value for a row
end

Base.@kwdef struct Seedset
    # type::Symbol
    filter::Vector{Term}
    change::Vector{Term}
    cnt::Int

    function Seedset(filter::Vector{Term}, change::Vector{Term}, cnt::Int) 
        allowed_filter_columns = [:cond, :status, :agegrp, :variant, :vaxstatus, :vaxrcvd, :tested, :quar]
        allowed_change_columns = [:cond, :status, :duration, :sickday, :variant, :vaxstatus, :vaxrcvd]

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


function Base.show(io::IO,  t::Term)
    if isa(t.val, Symbol)
        print(io, t.trait, "=:", t.val)
    else
        print(io, t.trait, "=", t.val)
    end
end

function Base.show(io::IO, v::Vector{Term})
    print(io, "[")
    for item in v
        show(io, item)
        print(io, ", ")
    end
    print(io, "]")
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
                if isempty(col[idx])
                    col[idx] = [ch.val]
                else
                    push!(col[idx], ch.val)
                end
            end
        else 
            col[filt] .= ch.val  # set the value for each trait column
        end
    end

    println("*** seed day ", DAY_CTR[:day], " count: ", s.cnt, " filter: ", s.filter, " change: ", s.change )

end


"""
    maketraitseedfunc( ; cond=:nil, variant=:base, duration=1, filter::Vector{Term}, cnt, forlocale=0, triggerdate, startofday)

Create a Seedset that contains the filter for whom to make sick, the traits to be set, and the cnt of people to be changed. 

**And** call seed\\_case\\_gen for you to return the callback function that encloses this Seedset.
"""
function maketraitseedfunc( ; cond=:nil, variant=:base, duration=1, filter::Vector{Term}, cnt, forlocale=0, triggerdate, forstartofday)
    ss = Seedset(filter=filter, cnt=cnt, 
            change=[Term(:status, :infectious), Term(:cond, cond), Term(:sickday, triggerdate), Term(:variant, variant), Term(:duration, duration)])

    seed_case_gen(ss; forlocale=forlocale, triggerdate=triggerdate, forstartofday=forstartofday)
end


"""
    makenotsickseedfunc( ; status, filter::Vector{Term}, cnt, forlocale=0, triggerdate, startofday)

Create a Seedset that contains the filter for whom to change from sick to either recovered or dead, and the cnt of people to be changed. 

This is just a shortcut for maketraitseedfunc with defaults to make people "well" again.

**And** call seed\\_case\\_gen for you to return the callback function that encloses this Seedset.
"""
function makenotsickseedfunc( ; status, filter::Vector{Term}, cnt, forlocale=0, triggerdate, forstartofday)
    if !in(status, [recovered, dead])
        throw(DomainError(status, "Status must be either recovered or dead"))
    end
    
    ss = Seedset(filter=filter, cnt=cnt,
            change=[Term(:status, status), Term(:cond, uninfected)])

    seed_case_gen(ss; forlocale=forlocale, triggerdate=triggerdate, forstartofday=forstartofday)
end

####################################################################
# functions to generate the callback functions for cases, specifically
#   seeding cases or setting sets of traits in the population
####################################################################

"""
Output a callback function that seeds the population (of a locale) with specific trait values. 
Sets a filter for which people are to be seeded with traits.
Sets the number of people to be seeded.
The callback function does NOT create people in the population. It changes the traits of people already in the population.

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(ss::Seedset; forlocale=0, triggerdate, forstartofday) # these args go into the called/returned seed! case

    function caserunner(locdat, socialparams, infectset, sdcases, age_idx_loc; day, startofday, locale)  # args must match runcases loop in run_a_sim
        if (day == triggerdate) & (startofday == forstartofday) 
            if (forlocale == 0) | (forlocale == locale)
                settraits!(locdat, ss::Seedset)
            end
        end
    end

end

"""
Generate seeding cases.
inputs: day, cnt, duration, cond, agegrp
Two of the inputs may refer to multiple items and must match in number of items.

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen_old(cnt, duration, cond, variant, agegrp; forlocale, triggerdate, forstartofday) # these args go into the returned seed! case
    # caserunner gets returned; assign it a value at the cmdline; use as an input to run_a_sim

    function caserunner(locdat, socialparams, infectset, sdcases, age_idx_loc; startofday, locale, day)  # args must match runcases loop in run_a_sim
        if (day == triggerdate) & (startofday == forstartofday) 
            if (forlocale == 0) | (forlocale == locale)
                seed!(cnt, duration, cond, variant, agegrp, locdat)  # payload: this is what the function will do when run
            end
        end
    end

end