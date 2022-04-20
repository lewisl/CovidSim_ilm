######################################################################################
# setup and initialization functions: ILM Model
######################################################################################


function setup(ndays, locales;  # must provide following inputs
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
        series = build_series_table(locales, agegrp, ndays)

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

    model = (ndays=ndays, locales=locales, dat=datadict, series=series, geo=geodata, 
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


Base.@kwdef struct Series
    cum::Dict{Int, Matrix{Int}} # locale as Int, matrix of cum history columns
    new::Dict{Int, Matrix{Int}} # locale as Int, matrix of new (each day) history columns
    groups::Vector{Symbol}
    cols::OrderedDict{Symbol, UnitRange{Int64}}
end


function build_series_table(locales, agegrp, n_days)

    cols = [Symbol(col,"_", age) for col in seriesgroups for age in vcat(collect(string.(instances(agegrp))),"total")]
    colvals = [zeros(Int,n_days) for _ in 1:length(cols)]
    series = Dict(loc => (cum = Table(; zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...), 
                          new = Table(; zip(cols,[zeros(Int,n_days) for _ in 1:length(cols)])...))
             for loc in locales)

    return series
end

function build_series_oldway(locales, n_days, map2dict=Dict())
    if isempty(map2dict)
        map2dict = OrderedDict{Symbol, UnitRange{Int64}}(
            :unexposed=>1:6, :infectious=>7:12, :recovered=>13:18, :dead=>19:24,            # status
            :nil=>25:30, :mild=>31:36, :sick=>37:42, :severe=>43:48, :totinfected=>49:54,   # condition
            :Pfizer=>55:60, :Moderna=>61:66, :JnJ=>67:72, :totvaccinated=>73:78,            # vaccine
            :base=>79:84, :alpha=>85:90, :delta=>91:96, :omicron_ba1=>97:102, :omicron_ba2=>103:108                # variant
            )   
    end
    group = collect(keys(map2dict))  

    tmpdict = Dict{Int, Series}()   # dict of locales to Series

    series = Series(groups = group,
                    cols = map2dict,
                    cum = Dict(loc => zeros(Int, n_days, map2dict[last(group)][end]) for loc in locales),
                    new = Dict(loc => zeros(Int, n_days, map2dict[last(group)][end]) for loc in locales)              
                )

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
    function build_spread_params(variantfilename, paramdir)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant.

Method to build infect params from one file that contains all variants.
"""
function build_infect_params(variantfilename, paramdir)
    infectdict = YAML.load_file(joinpath(paramdir, variantfilename), dicttype=Dict{Symbol, Any})

    infectset = build_spread_params(infectdict)
    (transitionset, trvec) = build_transition_params(infectdict)

    return infectset, transitionset, trvec
end


"""
    function build_spread_params(variants, paramdir)

Build parameters for the spread of infection and the immunity conferred by recovering
from infection for each variant.

Method to build spread params from a separate file for each variant.
"""
function build_spread_params(variants, paramdir)
    infectset = Dict{Symbol, Infectparams}()
    for variant in keys(variants)
        v = YAML.load_file(joinpath(paramdir, "variant_parameters", variants[variant][:directory_name],
            variants[variant][:infect_fname]), dicttype=Dict{Symbol, Any})
        v = Infectparams(v)
        infectset[variant] = v  # access a param as infectset[:alpha].recvrisk
    end

    return infectset
end


"""
    function build_spread_params(infectdict)

Method to build spread params from dict containing params for all variants.

This is the method model building actually uses!
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
        # if isempty(infectset[variant].sendrisk)
        #     append!(infectset[variant].sendrisk, infectset[:base].sendrisk .* infectset[variant].basemultiplier)
        # end
    end

    return infectset
end


"""
    function build_transition_params(variants, paramdir)

Build transition matrix from each illness condition to outcomes at each transition day
for someone who is infected.

Method for loading from transition params from a separate yaml file per each variant.

Returns (transitionset, trvec)
"""
function build_transition_params(variants, paramdir)
    transitionset = Dict()

    for variant in keys(variants)
        transitionset[variant] = setup_dt(joinpath(paramdir, "variant_parameters", variants[variant][:directory_name], 
            variants[variant][:transition_fname])) 
    end

    # pre-allocate trvec used in hot loop: no. of columns in transition array
    sz = size(first(first(transitionset[:base])[end])[end][:transition], 2)
    trvec = zeros(sz)

    return (transitionset, trvec)
end


"""
    function build_transition_params(infectdict)

This method loads all transition params for all variants from one dict, which contains
all variants.

This is the method model building actually uses!

Returns (transitionset, trvec)
"""
function build_transition_params(infectdict)
    loadvariants = keys(infectdict) # array of strings to array of symbols
    transitionset = Dict()

    for variant in loadvariants
        
        transitionset[Symbol(variant)] = Transitionparams(
            tree=(isnothing(infectdict[variant][:transition][:tree]) ? nothing : 
                    setup_dt(infectdict[variant][:transition][:tree])),
            factors=Transitionfactors(infectdict[variant][:transition][:factors])
            )

    end

    # pre-allocate trvec used in hot loop: no. of columns in transition array
    # sz = size(first(first(transitionset[:base])[end])[end][:transition], 2)
    sz = size(transitionset[:base].tree.age0_19[1].transition, 2)
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



