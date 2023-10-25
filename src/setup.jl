######################################################################################
# setup and initialization functions: ILM Model
######################################################################################

"""
Setup a model
Provides a definition of a model that can be saved and also allows re-running the simulation.

Pre-allocates all data storage for a simulation:
    day1: first calendar day
    ndays: number of days to run simulation
    geodata: data for each locale
    socialparams: parameters that affect transmission across people
    infectset: characteristics of each variant of the virus
    progressionset: parameters that affect how disease changes over time in an infected person
    variantlist: list of variants
    trvec: pre-allocated vector to hold on-the-fly calculated probabilities or progression to new disease condition
    vaxset: available vaccines and parameters for each vaccine
    vaxschedset: schedule for dispensing vaccines (net of vaccine resistant people)
    dat: a row for each person in a locale that tracks statistics for each person during the simulation
    series: "historical" statistics for outcomes at the end of each day(rows) of the simulation by new (change) and cumulative
    seriescolnames: column names for each statistic collected
"""
function setup(ndays::Int64, locales;  
    # must provide following inputs
    day1,
    dovax=false,
    paramdir,
    geofilename, 
    socialfilename,
    vaccinefilename,
    scheddir,
    variantfilename)


    geodata = buildgeodata(geofilename, paramdir)

    socialparams = build_socialparams(socialfilename, paramdir)

    # variants, spread parameters, progression arrays
    infectset, progressionset, trvec, variantlist = build_infect_params(variantfilename, paramdir)

    # vaccines  TODO this is not the right approach: test if we have vax inputs instead. Maybe?
    if dovax
        vaxset, vaxlist = build_vaxset(vaccinefilename, paramdir)
        vaxschedset = build_vaxschedset(scheddir, paramdir)
    else
        vaxset, vaxlist = Dict(), Symbol[]  # even "empty" needs to be typed correctly--empty what?
        vaxschedset = Dict()  # nothing
    end

    # simulation data matrix: rows = persons, columns = traits
    dat = build_data(locales, geodata, ndays)

    # history series columns and history series
    colgroups = [:statuscols=>Symbol.(STATUSES), :condcols=>push!(Symbol.(INFECTIOUS_CASES), :totinfected), 
                :vaxcols=>push!(Symbol.(vaxlist), :totvaccinated), :variantcols=>variantlist]

    seriescolnames = make_col_names_dict(colgroups)
    series = build_series_table(locales, ndays, day1, seriescolnames)

    # days that get indoor_uplift per locale for all days of the simulation
    indoor_seq = build_indoor_seq(locales, ndays, geodata, series, socialparams.indoor_uplift)

    model = (ndays=ndays, day1=day1, locales=locales, dat=dat, series=series, geo=geodata, 
            progressionset=progressionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset, 
            social=socialparams, trvec=trvec, variantlist=variantlist, vaxlist = vaxlist, 
            indoor_seq=indoor_seq, seriescolnames=seriescolnames)  

    return model
end

"""
    setup(yaml_model)

Create a complete simulation model from a previously saved YAML model definition. The YAML file must first be loaded with function yaml_to_model. 
The output model is identical to that created from input parameter files to the function buildsim. This output is a named tuple of all required model parameters. 
"""
function setup(yaml_model)
    ym = yaml_model

    day1 = Dates.Date(ym["day1"])
    dovax = ym["dovax"]
    ndays = ym["ndays"]
    locales = ym["locales"]

    geodata = buildgeodata(CSV.read(IOBuffer(ym["geofile"]), Table))

    dat = build_data(locales, geodata, ndays)

    socialparams = build_socialparams(YAML.load(ym["socialfile"], dicttype=OrderedDict{Symbol, Any}))

    # variants, spread parameters, progression arrays
    infectset, progressionset, trvec, variantlist = build_infect_params(YAML.load(ym["variantfile"], dicttype=Dict{Symbol, Any}))

    # vaccines  TODO this is not the right approach: test if we have vax inputs instead
    if dovax
        vaxset = build_vaxset(YAML.load(ym["vaccinefile"], dicttype=Dict{Symbol,Any}))
        vaxscheds = YAML.load(ym["vaxscheds"])  # a Dict{Any, Any}
        vaxschedset = build_vaxschedset(vaxscheds)
    else
        vaxset = Dict()  # nothing
        vaxschedset = Dict()  # nothing
    end

    series = build_series_table(ym["locales"], ym["ndays"], day1, seriescolnames) 

    # days that get indoor_uplift per locale for all days of the simulation
    indoor_seq = build_indoor_seq(locales, geodata, caldays, socialparams.indoor_uplift)

    model = (ndays=ndays, day1=day1, locales=locales, dat=dat, series=series, geo=geodata,
        progressionset=progressionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset,
        social=socialparams, trvec=trvec, variantlist=variantlist, vaxlist=vaxlist,
        indoor_seq=indoor_seq, seriescolnames=seriescolnames)

    return model
end


"""
Convert a vector of dates from a csv file in format "mm/dd/yyyy"
to a vector of Julia numeric Date values in format yyyy-mm-dd
"""
function quickdate(strdates)  # 20x faster than the built-in date parsing, e.g.--runs in 5% the time
    ret = [parse.(Int,i) for i in split.(strdates, '/')]
    ret = [Date.(i[3], i[1], i[2]) for i in ret]
end

"""
Pre-allocate and initialize population data for all locales in the simulation.
Calls pop_data for each locale.
"""
function build_data(locales, geodata, n_days)

    popdat = Dict(loc => pop_data(geodata.pop[geodata.fips .== loc][1]) for loc in locales)

    # precalculate agegrp indices = indices to rows for people in each age group
    agegrp_idx = Dict(loc => precalc_agegrp_filt(popdat[loc]).idx for loc in locales)
    
    return (popdat=popdat, agegrp_idx=agegrp_idx)
end


"""
    precalculate agegrp indices--these do not change during the simulation
"""
function precalc_agegrp_filt(dat)  # dat for a single locale
    agegrp_filt_bit = Dict(age => dat.agegrp .== age for age in AGEGRPS)
    agegrp_filt_idx = Dict(age => findall(agegrp_filt_bit[age]) for age in AGEGRPS)
    return (boolean=agegrp_filt_bit, idx=agegrp_filt_idx)
end


"""
Pre-allocate and initialize population data for one locale in the simulation.
Returns a TypedTable which is a tuple of arrays:
- each column is a trait of people
- each row is a person who lives in that locale
This table is updated each day of the simulation.
"""
function pop_data(pop; age_dist=AGE_DIST)
    
        parts = apportion(pop, age_dist)

        # must use comprehension to initialize vector of vector NOT fill--fill creates vectors at same address
        # LazyTable faster than Table; allows single rows to be modified
        dat = LazyTable(
            status = fill(:unexposed, pop),                                         
            agegrp = reduce(vcat,[fill(age, parts[i]) for (i,age) in enumerate(AGEGRPS)]), 
            cond = fill(:uninfected, pop),                                           
            duration = zeros(Int, pop),                                             
            variant = [Symbol[] for _ in 1:pop],                                   
            sickday = [Int[] for _ in 1:pop],                                                                   
            recovday = [Int[] for _ in 1:pop],                                 
            deadday = zeros(Int, pop),                                             
            ring = zeros(Int, pop),                                                
            sdcase = fill(:none, pop),                                          
            vaxstatus = fill(:none, pop),               # :none, :first, :full, :booster  maybe others later...
            vaxrcvd = [[:none] for _ in 1:pop],         # Vector{Symbol} of vaccine symbols  :Pfizer, :Moderna, :JnJ
            vaxday = [Int[] for _ in 1:pop],            
            tested = falses(pop),                                                   
            testday = zeros(Int, pop),                                              
            quar = falses(pop),                                                     
            quarday = zeros(Int, pop))                                             

    return dat       
end


"""
Pre-allocate and initialize table to hold history of the simulation, using LazyTables.
Returns a Dict of LazyTable with 2 keys:
- key cum is cumulative data for the entire locale. Or you may think of cum as the current value of a statistic.
- key new is the net change of a statistic for the entire locale. Note that this includes both additions and substractions. In other words, this
cannot be used as "new daily infections," for example.

For each table the structure is:
- columns are statistics that are recorded for each day of the simulation
- rows are days.

This table is updated at the end of each day of the simulation.
"""
function build_series_table(locales, n_days, day1, seriescolnames)
    caldays = range(day1, step=Day(1), length=n_days)

    cols = [col for group in values(seriescolnames) for item in values(group) for col in values(item)]
    colvals = [zeros(Int,n_days) for _ in 1:length(cols)]
    series = Dict(loc => (cum = LazyTable(; caldays=caldays, zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...), 
                          new = LazyTable(; caldays=caldays, zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...))
             for loc in locales)

    return series
end

# calculate which days get indoor_uplift for the entire simulation run instead of in a hot loop
function build_indoor_seq(locales, ndays, geodata, series, indoor_lift)
    indoor_seq = Dict(loc => ones(Float64, ndays) for loc in locales)

    for loc in locales
        caldays = series[loc].cum.caldays  # TODO: dumb because it's always the same, but difficul to unwrap

        indoor_end = Date(geodata.indoor_end[geodata.fips.==loc][1])
            year_end = year(indoor_end)
        indoor_start = Date(geodata.indoor_st[geodata.fips.==loc][1])
            year_start = year(indoor_start)

        if year_end == year_start  # start and end within a calendar year

            for i in eachindex(indoor_seq)
                testdate = Date(year_end, month(caldays[i]), day(caldays[i]))  # use relative year
                if (testdate >= indoor_start) & (testdate <= indoor_end) 
                    indoor_seq[loc][i] += indoor_lift
                end
            end

        elseif year_end > year_start  # start in first year, end in following year

            current_year = year(first(caldays))
            set_year = year_start

            for i in eachindex(indoor_seq[loc])
                if year(caldays[i]) > current_year
                    set_year = set_year == year_start ? year_end : year_start # toggle set_year
                    current_year = year(caldays[i])     # advance current_year
                end

                if (month(caldays[i]) == 2) & (day(caldays[i]) == 29)
                    continue  # the simulation year may be a leap year but the pseudo year is not
                end

                testdate = Date(set_year, month(caldays[i]), day(caldays[i]))
                if (testdate >= indoor_start) & (testdate <= indoor_end) 
                    indoor_seq[loc][i] *= indoor_lift
                end
            end

        else

            throw(DomainError((indoor_start_str, indoor_end_str), "Date for indoor_end must be > indoor_start"))

        end
    end
    return indoor_seq
end


# column names for series table returned as Dict
function make_col_names_dict(arr::Vector{Pair{Symbol, Vector{Symbol}}})     # Vector{Pair{Symbol, Vector}}
    agenames = collect((Symbol.(AGEGRPS)..., :total))
       #       group         item         age    colname
    ret = Dict{Symbol, Dict{Symbol, Dict{Symbol, Symbol}}}()
    for group in arr
        ret[group[1]] = gen_col_names_dict(group[2], agenames)
    end

    return ret
end


function gen_col_names_dict(items1, items2)
    Dict(zip(Symbol.(items1), [Dict(zip(items2, repeat_join([st], items2))) for st in items1]))
end


function buildgeodata(filename::String, paramdir)
    tmp = LazyTable(CSV.File(joinpath(paramdir, filename)))
    buildgeodata(tmp)
end

function buildgeodata(geotable::T) where T <: LazyTable
    LazyTable(geotable, 
        density_factor = shifter(geotable.density,0.9,1.25), 
        anchor         = quickdate(geotable.anchor),
        indoor_st      = quickdate(geotable.indoor_st),
        indoor_end     = quickdate(geotable.indoor_end)
        )
end




"""
    function build_infect_params(variantfilename, paramdir)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant. Build paramaters for progressioning infected people to
different conditions of the virus and to recover or die at the end.
"""
function build_infect_params(variantfilename, paramdir)

    variantdict = YAML.load_file(joinpath(paramdir, variantfilename), dicttype=Dict{Symbol, Any})

    build_infect_params(variantdict)
end


function build_infect_params(variantdict) 

    (infectset, variantlist) = build_spread_params(variantdict)
    (progressionset, trvec) = build_progression_params(variantdict)

    return infectset, progressionset, trvec, variantlist
end


"""
    function build_spread_params(variantdict)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant.
"""
function build_spread_params(variantdict::Dict)
    infectset = LittleDict{Symbol, Infectparams}()
    variantlist = collect(keys(variantdict))

    for variant in variantlist
        # result of merge only includes child keys of :spread and :immunity
        newdict = merge(variantdict[variant][:spread], variantdict[variant][:immunity])
        infectset[Symbol(variant)] = Infectparams(newdict)
    end


    # set recvrisk and sendrisk
    for variant in variantlist
        if variant === :base
            continue
        end

        # if no factors provided for this variant, apply multiplier to the base variant 
        if isempty(infectset[variant].recvrisk) & isempty(infectset[variant].sendrisk)  
            # use :base for both recvrisk and apply multiplier to sendrisk
            append!(infectset[variant].recvrisk, infectset[:base].recvrisk)
            append!(infectset[variant].sendrisk, infectset[:base].sendrisk .* infectset[variant].basemultiplier)
        elseif isempty(infectset[variant].recvrisk)         # use :base for recvrisk
            append!(infectset[variant].recvrisk, infectset[:base].recvrisk)
        else isempty(infectset[variant].sendrisk)           
            # use :base for sendrisk and apply multiplier
            append!(infectset[variant].sendrisk, infectset[:base].sendrisk .* infectset[variant].basemultiplier)
        end
    end

    return infectset, variantlist
end


"""
    function build_progression_params(variantdict)

This method loads all progression params for all variants from one dict, which contains
all variants.

Returns (progressionset, trvec)
"""
function build_progression_params(variantdict)
    variantlist = collect(keys(variantdict)) # array of strings to array of symbols
    progressionset = Dict{Symbol, ProgressionParams}()

    @assert :base in variantlist "Variants parameter file must contain a variant called :base--not there!"

    # build the progressionset for :base-->needed to build for other variants
    variant = :base
    @assert !isnothing(variantdict[variant][:progression_tree]) "progression tree for variant must be provided in parameter file--not there!"
    progressionset[Symbol(variant)] = ProgressionParams(
                                            tree=setup_dt(variantdict[variant][:progression_tree]),
                                            factors=ProgressionFactors(variantdict[variant][:progression_factors])
                                            )

    for variant in variantlist
        variant === :base && continue
        progressionset[Symbol(variant)] = ProgressionParams(
                tree=(  !isnothing(variantdict[variant][:progression_tree])   ?   
                            setup_dt(variantdict[variant][:progression_tree]) :    # progression tree was provided for this variant
                            setup_dt(deepcopy(progressionset[:base].tree), variantdict[variant][:progression_factors][:riskadjust])  # build the tree by adjusting :base
                     ),          
                factors=ProgressionFactors(variantdict[variant][:progression_factors]))
    end

    # pre-allocate trvec used in hot loop: no. of columns in progression array
    sz = size(progressionset[:base].tree.age0_19[5], 2)
    trvec = zeros(sz)
 
    return (progressionset, trvec)
end


function build_socialparams(socialfilename, paramdir)  # first step: read the input file

    social_inputs = YAML.load_file(joinpath(paramdir, socialfilename), dicttype=OrderedDict{Symbol, Any})

    build_socialparams(social_inputs)

end



function build_socialparams(social_inputs::T) where T <: AbstractDict  # build the data structures based on the inputs

    # check for all required params
        required_params = [:contactfactors, :touchfactors, :gammashape, :indoor_uplift]
        has_all = true
        lacking = []
        for p in required_params
            if !haskey(social_inputs, p)
                push!(lacking, p)
                has_all = false
            end
        end
        @assert has_all "required keys: $lacking missing"

        # build arrays for contactfactors and touchfactors
            # keys are agegrps
            # values are a dict of conditions with values = probabilities
        cfarr = zeros(length(keys(first(values(social_inputs[:contactfactors])))), length(keys(social_inputs[:contactfactors])))
        tfarr = zeros(length(keys(first(values(social_inputs[:touchfactors])))), length(keys(social_inputs[:touchfactors])))

        for (i, v1) in enumerate(sort(social_inputs[:contactfactors]))
            cfarr[:, i] .= Float64.(values(v1[2]))
        end
        for (i, v1) in enumerate(sort(social_inputs[:touchfactors]))
            tfarr[:, i] .= Float64.(values(v1[2]))
        end

    
    SocialParams(       # struct defined in CovidSim_ilm.jl
        gammashape      = Float64(social_inputs[:gammashape]),
        indoor_uplift   = Float64(social_inputs[:indoor_uplift]),
        contactfactors  = cfarr,
        touchfactors    = tfarr
        )
    
end


#####################################################################################
# Other helper functions
#####################################################################################

noop(args...; kwargs...) = nothing

function makemaptup(keys, values)
    NamedTuple{keys}(values)
end


#####################################################################################
# simple math helper functions
#   
#####################################################################################
"""
shifter makes linear changes in value ranges, preserving relative values
"""
@inline function shifter(x::AbstractArray, newmin, newmax)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
end

@inline @fastmath function shifter(x::Array, oldmin, oldmax, newmin, newmax)
    newmin .+ (newmax - newmin) / (oldmax - oldmin) .* (x .- oldmin)
end

@inline @fastmath function shifter(x::Float64, oldmin, oldmax, newmin, newmax)
    newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
end


@inline function shifter(x::AbstractArray, newval, mode::Symbol)
    if mode === :min
        newmin = newval
        newmax = maximum(x)
    elseif mode === :max
        newmin = minimum(x)
        newmax = newval
    else
        throw(DomainError(mode, "mode must be one of :min or :max"))
    end
    shifter(x, minimum(x), maximum(x), newmin, newmax)
end


@inline function shifter(x::Array; minmult=1.0, maxmult=1.0, mult=1.0)
    if mult != 1.0
        maxmult = minmult = mult
    end
    oldmin = minimum(x)
    oldmax = maximum(x)
    newmin = minmult * oldmin
    newmax = maxmult * oldmax
    shifter(x, oldmin, oldmax, newmin, newmax)
end


"""
    limdict(dct::Dict, op::Function)

Finds minimum or maxium value of the leaves of a dict.
Warning: not general! works on dict with 2 levels and 
numerical values at the lower level.
"""
function limdict(dct::AbstractDict, op::Function)
    minop = <
    cv = op == minop ? Inf : -Inf
    for v1 in values(dct)
        for v2 in values(v1)
            cv = op(v2, cv) ? v2 : cv
        end
    end
    return cv
end


"""
Warning: not general! works on dict with 2 levels and 
numerical values at the lower level.
"""
@inline function shifter(d::AbstractDict, newmin, newmax)
    ret = deepcopy(d)
    oldmin = limdict(d, <)
    oldmax = limdict(d, >)
    
    for k1 in keys(ret)
        for k2 in keys(ret[k1])
            x = ret[k1][k2]
            # ret[k1][k2] = newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
            ret[k1][k2] = shifter(x, oldmin, oldmax, newmin, newmax)
        end
    end

    return ret
end


function apportion(x::Int, splits::Array)  # x is the number to be split into portions
    @assert isapprox(sum(splits), 1.0)
    maxidx = argmax(splits)
    parts = round.(Int, splits .* x)
    diff = sum(parts) - x
    parts[maxidx] -= diff
    return parts
end
