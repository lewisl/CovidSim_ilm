#########################################################################################
# simstats.jl:    Simulation Stats -- very preliminary
#########################################################################################


gt(dict, key) = get(dict, key, 0)


function init_runtime_trace(locales, ndays)
    Dict(loc => (
            contacts = zeros(Int, ndays),
            touched = zeros(Int, ndays),
            spread_new_infected = zeros(Int, ndays),
            recovered = zeros(Int, ndays),
            dead = zeros(Int, ndays),
        ) for loc in locales)
end


function count_event_days!(dest, days)
    for day in days
        if 1 <= day <= length(dest)
            dest[day] += 1
        end
    end
    return dest
end


function count_event_days!(dest, day::Int)
    if 1 <= day <= length(dest)
        dest[day] += 1
    end
    return dest
end


function build_daily_trace(popdat, series, locale; runtime_trace=nothing, trace_ages=[:age60_79, :age80_up])
    locdat = popdat[locale]
    caldays = collect(series[locale].cum.caldays)
    ndays = length(caldays)

    contacts = isnothing(runtime_trace) ? zeros(Int, ndays) : copy(runtime_trace[locale].contacts)
    touched = isnothing(runtime_trace) ? zeros(Int, ndays) : copy(runtime_trace[locale].touched)
    spread_new_infected = isnothing(runtime_trace) ? zeros(Int, ndays) : copy(runtime_trace[locale].spread_new_infected)
    total_new_infected = zeros(Int, ndays)
    total_recovered = isnothing(runtime_trace) ? zeros(Int, ndays) : copy(runtime_trace[locale].recovered)
    total_dead = isnothing(runtime_trace) ? zeros(Int, ndays) : copy(runtime_trace[locale].dead)

    count_event_days!.(Ref(total_new_infected), locdat.sickday)
    if isnothing(runtime_trace)
        count_event_days!.(Ref(total_recovered), locdat.recovday)
        count_event_days!.(Ref(total_dead), locdat.deadday)
    end

    cols = Pair{Symbol, Any}[
        :day => collect(1:ndays),
        :calday => caldays,
        :contacts => contacts,
        :touched => touched,
        :new_infected => total_new_infected,
        :spread_new_infected => spread_new_infected,
        :recovered => total_recovered,
        :dead => total_dead,
    ]

    for age in trace_ages
        mask = locdat.agegrp .== age
        age_infected = zeros(Int, ndays)
        age_recovered = zeros(Int, ndays)
        age_dead = zeros(Int, ndays)

        count_event_days!.(Ref(age_infected), locdat.sickday[mask])
        count_event_days!.(Ref(age_recovered), locdat.recovday[mask])
        count_event_days!.(Ref(age_dead), locdat.deadday[mask])

        push!(cols, Symbol(:new_infected_, age) => age_infected)
        push!(cols, Symbol(:recovered_, age) => age_recovered)
        push!(cols, Symbol(:dead_, age) => age_dead)
    end

    LazyTable(; (name => value for (name, value) in cols)...)
end


function runsim_trace(model; trace_ages=[:age60_79, :age80_up], kwargs...)
    runtime_trace = init_runtime_trace(model.locales, model.ndays)
    popdat, series = runsim(model; runtime_trace=runtime_trace, kwargs...)
    trace = Dict(loc => build_daily_trace(popdat, series, loc; runtime_trace=runtime_trace, trace_ages=trace_ages)
                 for loc in model.locales)
    return popdat, series, trace
end


function runsim_spread_debug(model; config=SpreadDebugConfig(), kwargs...)
    spread_debug_trace = init_spread_debug_trace(model.locales; config=config)
    popdat, series = runsim(model; spread_debug_trace=spread_debug_trace, kwargs...)
    trace = Dict(loc => spread_debug_tables(spread_debug_trace[loc]) for loc in model.locales)
    return popdat, series, trace
end


"""
    function stat1(series, locale)

Returns a table of status outcomes for the simulation by age group.
"""
function stat_cond(series, locale)
    cumhist = series[locale].cum
    newhist = series[locale].new

    stat1 = Table(item=Symbol[], total=Int[], age0_19=Int[], age20_39=Int[], age40_59=Int[], age60_79=Int[], age80_up=Int[])
    calc_cols = [:total, :age0_19, :age20_39, :age40_59, :age60_79, :age80_up]

    # TODO use list of ages rather than hardcoding
    stat1 = push!(stat1, (item=:pop,
                    total    = getproperty(cumhist, :unexposed_total)[1] + getproperty(cumhist, :infectious_total)[1],
                    age0_19  = getproperty(cumhist, :unexposed_age0_19)[1] + getproperty(cumhist, :infectious_age0_19)[1],
                    age20_39 = getproperty(cumhist, :unexposed_age20_39)[1] + getproperty(cumhist, :infectious_age20_39)[1],
                    age40_59 = getproperty(cumhist, :unexposed_age40_59)[1] + getproperty(cumhist, :infectious_age40_59)[1],
                    age60_79 = getproperty(cumhist, :unexposed_age60_79)[1] + getproperty(cumhist, :infectious_age60_79)[1],
                    age80_up = getproperty(cumhist, :unexposed_age80_up)[1] + getproperty(cumhist, :infectious_age80_up)[1])
                )
                
    push!(stat1, (item=:died,
                  total    = getproperty(cumhist, :dead_total)[end],
                  age0_19  = getproperty(cumhist, :dead_age0_19)[end],
                  age20_39 = getproperty(cumhist, :dead_age20_39)[end],
                  age40_59 = getproperty(cumhist, :dead_age40_59)[end],
                  age60_79 = getproperty(cumhist, :dead_age60_79)[end],
                  age80_up = getproperty(cumhist, :dead_age80_up)[end])
                  )

    push!(stat1, (item=:unexposed,
                  total    = getproperty(cumhist, :unexposed_total)[end],
                  age0_19  = getproperty(cumhist, :unexposed_age0_19)[end],
                  age20_39 = getproperty(cumhist, :unexposed_age20_39)[end],
                  age40_59 = getproperty(cumhist, :unexposed_age40_59)[end],
                  age60_79 = getproperty(cumhist, :unexposed_age60_79)[end],
                  age80_up = getproperty(cumhist, :unexposed_age80_up)[end])
                  )

    # inter-row calculation infected = pop - unexposed
    poprow = stat1[stat1.item .=== :pop]
    unexprow = stat1[stat1.item .=== :unexposed]
    newrow = zeros(Int, 6)
    for (i, col) in enumerate(calc_cols)
        newrow[i] = getproperty(poprow, col)[] - getproperty(unexprow, col)[]
    end

    push!(stat1, (item = :ever_infected,
            total = newrow[1], age0_19 = newrow[2], age20_39 = newrow[3], 
            age40_59 = newrow[4], age60_79 = newrow[5], age80_up = newrow[6])
        )

    # inter-row calculation recovered = infected - died
    diedrow = stat1[stat1.item .=== :died]
    infectedrow = stat1[stat1.item .=== :ever_infected]
    newrow[:] .= 0
    for (i, col) in enumerate(calc_cols)
        newrow[i] = (getproperty(infectedrow, col)[] - getproperty(diedrow, col)[] 
                    - getproperty(cumhist, Symbol(:infectious, "_", col))[end])
    end

    push!(stat1, (item=:recovered,
                    total    = newrow[1], age0_19  = newrow[2], age20_39 = newrow[3],
                    age40_59 = newrow[4], age60_79 = newrow[5], age80_up = newrow[6])
                  )



    # add pct of population columns
    statpct = Table(item = stat1.item, total = stat1.total, total_pct = stat1.total ./ stat1.total[1],   # 
            age0_19 = stat1.age0_19, age0_19_pct = stat1.age0_19 ./ stat1.age0_19[1], 
            age20_39=stat1.age20_39, age20_39_pct = stat1.age20_39 ./ stat1.age20_39[1],
            age40_59 = stat1.age40_59, age40_59_pct = stat1.age40_59 ./ stat1.age40_59[1], 
            age60_79=stat1.age60_79, age60_79_pct = stat1.age60_79 ./ stat1.age60_79[1], 
            age80_up=stat1.age80_up, age80_up_pct = stat1.age80_up ./ stat1.age80_up[1])

    return statpct

end


function stat_vax(popdat, locale)
    td = @Select(vaxrcvd, agegrp, status)(popdat[locale]) # view with 2 columns on popdat
    vaxes = td.vaxrcvd

    stat_vax = Table(item=Symbol[], total=Int[], age0_19=Int[], age20_39=Int[], age40_59=Int[], age60_79=Int[], age80_up=Int[])
    calc_cols = [:total, :age0_19, :age20_39, :age40_59, :age60_79, :age80_up]
    
    col_age0_19  = countmap(last.(vaxes[(td.agegrp .== :age0_19) .& (td.status .!= :dead)]))
    col_age20_39 = countmap(last.(vaxes[(td.agegrp .== :age20_39) .& (td.status .!= :dead)]))
    col_age40_59 = countmap(last.(vaxes[(td.agegrp .== :age40_59) .& (td.status .!= :dead)]))
    col_age60_79 = countmap(last.(vaxes[(td.agegrp .== :age40_59) .& (td.status .!= :dead)]))
    col_age80_up = countmap(last.(vaxes[(td.agegrp .== :age80_up) .& (td.status .!= :dead)]))
    col_total    = countmap(last.(vaxes))
    
    for k in (:none, :Moderna, :Pfizer, :JnJ)
        push!(stat_vax, 
              (item=k, total=gt(col_total,k), age0_19=gt(col_age0_19,k), age20_39=gt(col_age20_39,k), age40_59=gt(col_age40_59,k),
               age60_79=gt(col_age60_79,k), age80_up=gt(col_age80_up,k))
               )
    end

    stat_pct = @Select(item, total, total_pct = $total ./ sum($total),
                age0_19,  age0_19_pct  = $age0_19  ./ sum($age0_19),
                age20_39, age20_39_pct = $age20_39 ./ sum($age20_39),
                age40_59, age40_59_pct = $age40_59 ./ sum($age40_59),
                age60_79, age60_79_pct = $age60_79 ./ sum($age60_79),
                age80_up, age80_up_pct = $age80_up ./ sum($age80_up))(stat_vax)

    return stat_pct

end


function stat_repeat(popdat, locale)
    thisdat = @Select(variant, agegrp)(popdat[locale]) # view with 2 columns on popdat
    variants = thisdat.variant

    stat_count = Table(item=Symbol[], total=Int[], age0_19=Int[], age20_39=Int[], age40_59=Int[], age60_79=Int[], age80_up=Int[])

    col_total    = sort(countmap(length.(variants)))
    col_age0_19  = sort(countmap(length.(variants[thisdat.agegrp .== :age0_19])))
    col_age20_39 = sort(countmap(length.(variants[thisdat.agegrp .== :age20_39])))
    col_age40_59 = sort(countmap(length.(variants[thisdat.agegrp .== :age40_59])))
    col_age60_79 = sort(countmap(length.(variants[thisdat.agegrp .== :age40_59])))
    col_age80_up = sort(countmap(length.(variants[thisdat.agegrp .== :age80_up])))

    for k in 0:10
        push!(stat_count, 
              (item=Symbol("Times_",k), total=gt(col_total,k), age0_19=gt(col_age0_19,k), age20_39=gt(col_age20_39,k), age40_59=gt(col_age40_59,k),
               age60_79=gt(col_age60_79,k), age80_up=gt(col_age80_up,k))
               )
    end

    stat_pct = @Select(item, total, total_pct = $total ./ sum($total),
                        age0_19,  age0_19_pct  = $age0_19  ./ sum($age0_19),
                        age20_39, age20_39_pct = $age20_39 ./ sum($age20_39),
                        age40_59, age40_59_pct = $age40_59 ./ sum($age40_59),
                        age60_79, age60_79_pct = $age60_79 ./ sum($age60_79),
                        age80_up, age80_up_pct = $age80_up ./ sum($age80_up))(stat_count)

    return stat_pct
end


function stat_breakout(popdat, locale)

    @Select(vaxday, sickday, agegrp)(locdat[((locdat.status .== infectious) .| (locdat.status .== recovered)) .& 
    (locdat.vaxstatus .!= :none) .& (last.(locdat.sickday) .> first.(locdat.vaxday))])

end


function vax_summary(model, locale)
    cumhist = model.series[locale].cum
    newhist = model.series[locale].new
    locdat = model.dat[locale]
    day1 = model.day1
    n = length(newhist)

    # construct the data series: rows are vax sequence
    daily_vaxes = Table(zeros(Int, n, numcols))
    maxtimes = maximum(length.(locdat.sickday))  # maximum no. of times anyone has gotten infected
    for i in 2:maxtimes+1
        for age in AGEGRPS  # accumulate all days on which anyone got sick the 1st, 2nd, 3rd... time
            daygotsick = countmap(get.(locdat.sickday[locdat.agegrp .== age],i,0))
            for (k,v) in daygotsick
                if k == 0
                    continue  # ignore people who never got infected
                end
                daily_cases_series[k, Int(age)] += v
            end
        end
    end

end

function detailed_vax_series(model, locale)
    cumhist = model.series[locale].cum
    newhist = model.series[locale].new
    locdat = model.dat[locale]
    day1 = model.day1
    n = length(newhist)

    # construct the data series
    daily_vaxes = Table(zeros(Int, n, numcols))
    maxtimes = maximum(length.(locdat.sickday))  # maximum no. of times anyone has gotten infected
    for i in 2:maxtimes+1
        for age in AGEGRPS  # accumulate all days on which anyone got sick the 1st, 2nd, 3rd... time
            daygotsick = countmap(get.(locdat.sickday[locdat.agegrp .== age],i,0))
            for (k,v) in daygotsick
                if k == 0
                    continue  # ignore people who never got infected
                end
                daily_cases_series[k, Int(age)] += v
            end
        end
    end
end

function daily_cases_series!(model, locale)
    cumhist = model.series[locale].cum
    newhist = model.series[locale].new
    locdat = model.dat[locale]
    scn = seriescolnames = model.seriescolnames
    day1 = model.day1
    n = length(newhist)

    # construct the data series
    daily_cases = Table(zeros(Int, n, numcols))
    maxtimes = maximum(length.(locdat.sickday))  # maximum no. of times anyone has gotten infected
    for i in 2:maxtimes+1
        for age in AGEGRPS  # accumulate all days on which anyone got sick the 1st, 2nd, 3rd... time
            daygotsick = countmap(get.(locdat.sickday[locdat.agegrp .== age],i,0))
            for (k,v) in daygotsick
                if k == 0
                    continue  # ignore people who never got infected
                end
                daily_cases_series[k, Int(age)] += v
            end
        end
    end

    daily_cases_series[:, numcols] = sum(daily_cases_series, dims=2)  # sum the columns across each row
end
