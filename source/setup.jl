######################################################################################
# setup and initialization functions: ILM Model
######################################################################################


function setup(ndays, locales;  # must provide following inputs
    day1,
    dovax=false,
    dovariant=false,
    paramdir,
    geofilename, 
    socialfilename,
    vaccinefilename,
    variantfilename)

    # geodata
        geodata = buildgeodata(geofilename)

    # simulation data matrix
        datadict = build_data(locales, geodata, ndays)

    # history series
        series = build_series_table(locales, agegrp, ndays, day1)

    # social parameters
        socialparams = build_socialparams(socialfilename, paramdir)

    # variants, spread parameters, transition arrays
        infectset, transitionset, trvec = build_infect_params(variantfilename, paramdir)


    # vaccines  TODO this is not the right approach: test if we have vax inputs instead
    if dovax
        vaxset = build_vaxset(vaccinefilename, paramdir=paramdir)
        vaxschedset = build_vaxschedset()
    else
        vaxset = Dict()  # nothing
        vaxschedset = Dict()  # nothing
    end

    model = (ndays=ndays, day1=day1, locales=locales, dat=datadict, series=series, geo=geodata, 
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

    pop = [geodata[geodata[:, "fips"] .== loc, "pop"][1] for loc in locales]

    popdat = Dict(loc => pop_data(geodata[geodata[:, "fips"] .== loc, "pop"][1]) for loc in locales)

    # precalculate agegrp indices
    agegrp_idx = Dict(loc => precalc_agegrp_filt(popdat[loc]).idx for loc in locales)
    
    return Dict("popdat"=>popdat, "agegrp_idx"=>agegrp_idx)
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
            vaxstatus = fill(:none, pop),          # :none, :first, :multiple, :full, :booster  maybe others later...
            vaxrcvd = [[:none] for _ in 1:pop],    # Vector{Symbol} of vaccine symbols  :Pfizer, :Moderna, :JnJ
            vaxday = [[0] for _ in 1:pop],                                          # Vector{Int}
            tested = falses(pop),                                                   # Bool
            testday = zeros(Int, pop),                                              # Vector{Int}
            quar = falses(pop),                                                     # Bool
            quarday = zeros(Int, pop))                                              # Int

    return dat       
end


function build_series_table(locales, agegrp, n_days, day1)
    calday = range(day1, step=Day(1), length=n_days)
    cols = [Symbol(col,"_", age) for col in seriesgroups for age in vcat(collect(string.(instances(agegrp))),"total")]
    colvals = [zeros(Int,n_days) for _ in 1:length(cols)]
    series = Dict(loc => (cum = Table(; calday=calday, zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...), 
                          new = Table(; calday=calday, zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...))
             for loc in locales)

    return series
end


function buildgeodata(filename)
    geo = DataFrame(CSV.File(filename))
    insertcols!(geo, "density_factor" => shifter(geo[:, "density"],0.9,1.25))

    # fix dates   
    insertcols!(geo, "anchor2" => quickdate(geo[:, "anchor"]))
    insertcols!(geo, "limit2" => quickdate(geo[:, "limit"]))
    select!(geo, Not(["anchor", "limit"]))
    rename!(geo, "anchor2" => "anchor")
    rename!(geo, "limit2" => "limit")

    return geo
end

"""
    function build_infect_params(variantfilename, paramdir)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant. Build paramaters for transitioning infected people to
different conditions of the virus and to recover or die at the end.
"""
function build_infect_params(variantfilename, paramdir)
    infectdict = YAML.load_file(joinpath(paramdir, variantfilename), dicttype=Dict{Symbol, Any})

    infectset = build_spread_params(infectdict)
    (transitionset, trvec) = build_transition_params(infectdict)

    return infectset, transitionset, trvec
end


"""
    function build_spread_params(infectdict)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant.
"""
function build_spread_params(infectdict::Dict)
    infectset = LittleDict{Symbol, Infectparams}()
    loadvariants = collect(keys(infectdict))
    if isempty(variantlist)
        append!(variantlist, loadvariants) 
    end

    for variant in loadvariants
        newdict = merge(infectdict[variant][:spread], infectdict[variant][:immunity])
        # newdict = Dict(Symbol(k) => v for (k,v) in newdict)
        infectset[Symbol(variant)] = Infectparams(newdict)
    end


    # set recvrisk and sendrisk
    for variant in loadvariants
        if variant == :base
            continue
        end
        if isempty(infectset[variant].recvrisk)
            append!(infectset[variant].recvrisk, infectset[:base].recvrisk .* infectset[variant].basemultiplier)
            append!(infectset[variant].sendrisk, infectset[:base].sendrisk)
        end
    end

    return infectset
end


"""
    function build_transition_params(infectdict)

This method loads all transition params for all variants from one dict, which contains
all variants.

Returns (transitionset, trvec)
"""
function build_transition_params(infectdict)
    loadvariants = keys(infectdict) # array of strings to array of symbols
    transitionset = Dict{Symbol, Transitionparams}()

    @assert :base in loadvariants "Variants parameter file must contain a variant called :base--not there!"

    # build the transitionset for :base-->needed to build for other variants
    variant = :base
    @assert !isnothing(infectdict[variant][:transition][:tree]) "transition tree for variant must be provided in parameter file--not there!"
    transitionset[Symbol(variant)] = Transitionparams(
                                            tree=setup_dt(infectdict[variant][:transition][:tree]),
                                            factors=Transitionfactors(infectdict[variant][:transition][:factors])
                                            )

    for variant in loadvariants
        variant == :base && continue
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
        @assert has_all "required keys: $lacking not in $(infectfilename)"

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
