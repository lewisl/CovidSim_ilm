######################################################################################
# setup and initialization functions: ILM Model
######################################################################################


function setup(ndays, locales;  # must provide following inputs
    day1,
    dovax=false,
    paramdir,
    geofilename, 
    socialfilename,
    vaccinefilename,
    scheddir,
    variantfilename)

    # geodata
        geodata = buildgeodata(geofilename)

    # social parameters
        socialparams = build_socialparams(socialfilename, paramdir)

    # variants, spread parameters, transition arrays
        infectset, transitionset, trvec, variantlist = build_infect_params(variantfilename, paramdir)

    # vaccines  TODO this is not the right approach: test if we have vax inputs instead
    if dovax
        vaxset, vaxlist = build_vaxset(vaccinefilename, paramdir)
        vaxschedset = build_vaxschedset(scheddir, paramdir)
    else
        vaxset, vaxlist = Dict(), []  # nothing
        vaxschedset = Dict()  # nothing
    end

    # simulation data matrix
    dat = build_data(locales, geodata, ndays)

    # history series columns and history series
        colgroups = [:statuscols=>statuses, :condcols=>push!(Symbol.(infectious_cases), :totinfected), 
                    :vaxcols=>push!(Symbol.(vaxlist), :totvaccinated), :variantcols=>variantlist]
        seriescolnames = make_col_names_dict(colgroups)
        series = build_series_table(locales, agegrp, ndays, day1, seriescolnames)

    model = (ndays=ndays, day1=day1, locales=locales, dat=dat, series=series, geo=geodata, 
            transitionset=transitionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset, 
            social=socialparams, trvec=trvec, variantlist=variantlist, vaxlist = vaxlist, 
            seriescolnames=seriescolnames)  

    return model
end

"""
    setup(yaml_model)

Create a complete simulation model from a previously saved YAML model definition. The YAML file must first be loaded with function yaml_to_model. The output model is identical to that created from input parameter files to the function buildsim. This output is a named tuple of all required model parameters. 
"""
function setup(yaml_model)
    ym = yaml_model

    day1 = Dates.Date(ym["day1"])
    dovax = ym["dovax"]
    ndays = ym["ndays"]
    locales = ym["locales"]

    #geodata
        geodata = buildgeodata(CSV.read(IOBuffer(ym["geofile"]), Table))

    # simulation data matrix
        dat = build_data(locales, geodata, ndays)

        
    # social parameters
        socialparams = build_socialparams(YAML.load(ym["socialfile"], dicttype=OrderedDict{Symbol, Any}))

    # variants, spread parameters, transition arrays
        infectset, transitionset, trvec, variantlist = build_infect_params(YAML.load(ym["variantfile"], dicttype=Dict{Symbol, Any}))


    # vaccines  TODO this is not the right approach: test if we have vax inputs instead
    if dovax
        vaxset = build_vaxset(YAML.load(ym["vaccinefile"], dicttype=Dict{Symbol,Any}))
        vaxscheds = YAML.load(ym["vaxscheds"])  # a Dict{Any, Any}
        vaxschedset = build_vaxschedset(vaxscheds)
    else
        vaxset = Dict()  # nothing
        vaxschedset = Dict()  # nothing
    end

    # history series
    series = build_series_table(ym["locales"], agegrp, ym["ndays"], day1) 

    model = (ndays=ndays, day1=day1, locales=locales, dat=dat, series=series, geo=geodata, 
            transitionset=transitionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset, 
            social=socialparams, trvec=trvec)  

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

    # pop = [geodata[geodata[:, "fips"] .== loc, "pop"][1] for loc in locales]

    popdat = Dict(loc => pop_data(geodata.pop[geodata.fips .== loc][1]) for loc in locales)

    # precalculate agegrp indices
    agegrp_idx = Dict(loc => precalc_agegrp_filt(popdat[loc]).idx for loc in locales)
    
    return (popdat=popdat, agegrp_idx=agegrp_idx)
end


"""
    precalculate agegrp indices--these do not change during the simulation
"""
function precalc_agegrp_filt(dat)  # dat for a single locale
    agegrp_filt_bit = Dict(age => dat.agegrp .== age for age in agegrps)
    agegrp_filt_idx = Dict(age => findall(agegrp_filt_bit[age]) for age in agegrps)
    return (boolean=agegrp_filt_bit, idx=agegrp_filt_idx)
end


"""
Pre-allocate and initialize population data for one locale in the simulation.
Returns a TypedTable which is a tuple of arrays:
- each column is a trait of people
- rows are days of the simulatoin
"""
function pop_data(pop; age_dist=age_dist)
        parts = apportion(pop, age_dist)

        # must use comprehension to initialize vector of vector NOT fill--fill creates identical vectors
        dat = Table(
            status = fill(unexposed, pop),                                          # enum status
            agegrp = reduce(vcat,[fill(age, parts[Int(age)]) for age in agegrps]),  # enum agegrp
            cond = fill(uninfected, pop),                                           # enum condition
            duration = zeros(Int, pop),                                             # Int
            variant = [Symbol[] for _ in 1:pop],                                    # Vector{Symbol}
            sickday = [[0] for _ in 1:pop],                                         # Vector{Int}
            recovday = [[0] for _ in 1:pop],                                        # Vector{Int}
            deadday = zeros(Int, pop),                                              # Int
            ring = zeros(Int, pop),                                                 # Int (not used as yet)
            sdcomply = fill(:none, pop),                                            # Symbol
            vaxstatus = fill(:none, pop),          # :none, :first, :full, :booster  maybe others later...
            vaxrcvd = [[:none] for _ in 1:pop],    # Vector{Symbol} of vaccine symbols  :Pfizer, :Moderna, :JnJ
            vaxday = [[0] for _ in 1:pop],                                          # Vector{Int}
            tested = falses(pop),                                                   # Bool
            testday = zeros(Int, pop),                                              # Vector{Int}
            quar = falses(pop),                                                     # Bool
            quarday = zeros(Int, pop))                                              # Int

    return dat       
end


function build_series_table(locales, agegrp, n_days, day1, seriescolnames)
    calday = range(day1, step=Day(1), length=n_days)
    # cols = [col for group in seriescolnames for item in group for col in item]
    cols = [col for group in values(seriescolnames) for item in values(group) for col in values(item)]
    colvals = [zeros(Int,n_days) for _ in 1:length(cols)]
    series = Dict(loc => (cum = Table(; calday=calday, zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...), 
                          new = Table(; calday=calday, zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...))
             for loc in locales)

    return series
end


# column names for series table returned as Dict
function make_col_names_dict(arr::Vector{Pair{Symbol, Vector}})
    agenames = collect((Symbol.(agegrps)..., :total))
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


function buildgeodata(filename::String)
    tmp = Table(CSV.File(filename))
    buildgeodata(tmp)
end

function buildgeodata(geotable::T) where T <: Table
    Table(geotable, 
        density_factor = shifter(geotable.density,0.9,1.25), 
        anchor         = quickdate(geotable.anchor),
        limit          = quickdate(geotable.limit)
        )
end




"""
    function build_infect_params(variantfilename, paramdir)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant. Build paramaters for transitioning infected people to
different conditions of the virus and to recover or die at the end.
"""
function build_infect_params(variantfilename, paramdir)

    infectdict = YAML.load_file(joinpath(paramdir, variantfilename), dicttype=Dict{Symbol, Any})

    build_infect_params(infectdict)
end


function build_infect_params(infectdict) 

    (infectset, variantlist) = build_spread_params(infectdict)
    (transitionset, trvec) = build_transition_params(infectdict)

    return infectset, transitionset, trvec, variantlist
end


"""
    function build_spread_params(infectdict)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant.
"""
function build_spread_params(infectdict::Dict)
    infectset = LittleDict{Symbol, Infectparams}()
    variantlist = collect(keys(infectdict))

    for variant in variantlist
        newdict = merge(infectdict[variant][:spread], infectdict[variant][:immunity])
        infectset[Symbol(variant)] = Infectparams(newdict)
    end


    # set recvrisk and sendrisk
    for variant in variantlist
        if variant === :base
            continue
        end
        if isempty(infectset[variant].recvrisk) & isempty(infectset[variant].sendrisk)  # use :base for both recvrisk and sendrisk
            append!(infectset[variant].recvrisk, infectset[:base].recvrisk .* infectset[variant].basemultiplier)
            append!(infectset[variant].sendrisk, infectset[:base].sendrisk)
        elseif isempty(infectset[variant].recvrisk)         # use :base for recvrisk
            append!(infectset[variant].recvrisk, infectset[:base].recvrisk .* infectset[variant].basemultiplier)
        else isempty(infectset[variant].sendrisk)           # use :base for sendrisk
            append!(infectset[variant].sendrisk, infectset[:base].sendrisk .* infectset[variant].basemultiplier)
        end
    end

    return infectset, variantlist
end


"""
    function build_transition_params(infectdict)

This method loads all transition params for all variants from one dict, which contains
all variants.

Returns (transitionset, trvec)
"""
function build_transition_params(infectdict)
    variantlist = collect(keys(infectdict)) # array of strings to array of symbols
    transitionset = Dict{Symbol, Transitionparams}()

    @assert :base in variantlist "Variants parameter file must contain a variant called :base--not there!"

    # build the transitionset for :base-->needed to build for other variants
    variant = :base
    @assert !isnothing(infectdict[variant][:transition][:tree]) "transition tree for variant must be provided in parameter file--not there!"
    transitionset[Symbol(variant)] = Transitionparams(
                                            tree=setup_dt(infectdict[variant][:transition][:tree]),
                                            factors=Transitionfactors(infectdict[variant][:transition][:factors])
                                            )

    for variant in variantlist
        variant === :base && continue
        transitionset[Symbol(variant)] = Transitionparams(
                tree=(  !isnothing(infectdict[variant][:transition][:tree])   ?   
                            setup_dt(infectdict[variant][:transition][:tree]) :    # transition tree was provided for this variant
                            setup_dt(deepcopy(transitionset[:base].tree), infectdict[variant][:transition][:factors][:riskadjust])  # build the tree by adjusting :base
                     ),          
                factors=Transitionfactors(infectdict[variant][:transition][:factors]))
    end

    # pre-allocate trvec used in hot loop: no. of columns in transition array
    sz = size(transitionset[:base].tree.age0_19[5], 2)
    trvec = zeros(sz)
 
    return (transitionset, trvec)
end


function build_socialparams(socialfilename, paramdir)

    social_inputs = YAML.load_file(joinpath(paramdir, socialfilename), dicttype=OrderedDict{Symbol, Any})

    build_socialparams(social_inputs)

end



function build_socialparams(social_inputs::T) where T <: AbstractDict

    # social_inputs = YAML.load_file(joinpath(paramdir, socialfilename), dicttype=OrderedDict{Symbol, Any})

    # check for all required params
        required_params = [:contactfactors, :touchfactors, :gammashape]
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


    Socialparams(
        gammashape      = Float64(social_inputs[:gammashape]),
        contactfactors  = cfarr,
        touchfactors    = tfarr
        )
    
end


#####################################################################################
# Other helper functions
#####################################################################################

noop(args...; kwargs...) = nothing


function make_an_enum!(name, strarr; pr=false)
    eval(:(@enum $(Symbol(name)) $(Symbol.(strarr)...)))  
    if pr
        
        display("text/markdown",  """**Created this enum** \n
        """)
        eval(Symbol(name)) 
        
    end
end

function make_an_enum!(name, strarr, start; pr=false)  # method with start value
    eval(:(@enum $(Symbol(name)) ($(Symbol(strarr[1])) = $(start)) $(Symbol.(strarr[2:end])...)  ))   

    if pr
        
        display("text/markdown",  """**Created this enum** \n
        """)
        eval(Symbol(name)) 
        
    end
end

function makemaptup(keys, values)
    NamedTuple{keys}(values)
end


#####################################################################################
# dodgy math helper functions
#####################################################################################

@inline @fastmath function shifter(x::Array, oldmin, oldmax, newmin, newmax)
    newmin .+ (newmax - newmin) / (oldmax - oldmin) .* (x .- oldmin)
end

@inline @fastmath function shifter(x::Float64, oldmin, oldmax, newmin, newmax)
    newmin + (newmax - newmin) / (oldmax - oldmin) * (x - oldmin)
end

@inline function shifter(x::Array, newmin, newmax)
    oldmin = minimum(x)
    oldmax = maximum(x)
    shifter(x, oldmin, oldmax, newmin, newmax)
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


function apportion(x::Int, splits::Array)
    @assert isapprox(sum(splits), 1.0)
    maxidx = argmax(splits)
    parts = round.(Int, splits .* x)
    diff = sum(parts) - x
    parts[maxidx] -= diff
    return parts
end
