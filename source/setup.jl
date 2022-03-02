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
        series = build_series(locales, ndays)

    # social parameters
        socialparams = build_socialparams(socialfilename, paramdir)

    # variants, spread parameters, transition arrays
        infectset, transitionset, trvec = build_infect_params(variantfilename, paramdir)


    # vaccines  TODO this is not the right approach: test if we have vax inputs instead
    if dovax
        vaxset = build_vaxset(vaccinefilename, paramdir=paramdir)
        vaxschedset = build_vaxschedset()
    else
        vaxset = nothing
        vaxschedset = nothing
    end

    return (ndays=ndays, locales=locales, dat=datadict, series=series, geo=geodata, 
            transitionset=transitionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset, 
            social=socialparams, trvec=trvec)  
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


# columns of history series: traits by agegrp and total:  first 5 cols are agegrps, 6th is total
# const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24,          # status
#                     nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54,   # conditions
#                     Pfizer=55:60, Moderna=61:66, JnJ=67:72, totvaccinated=73:78,          # vaccines
#                     base=79:84, alpha=85:90, delta=91:96, omicron=97:102)                 # variants


Base.@kwdef struct Series
    cum::Matrix{Int}
    new::Matrix{Int}
    groups::Vector{Symbol}
    cols::OrderedDict{Symbol, UnitRange{Int64}}
end


function build_series(locales, n_days)
    map2dict = OrderedDict{Symbol, UnitRange{Int64}}(
        :unexposed=>1:6, :infectious=>7:12, :recovered=>13:18, :dead=>19:24,            # status
        :nil=>25:30, :mild=>31:36, :sick=>37:42, :severe=>43:48, :totinfected=>49:54,   # condition
        :Pfizer=>55:60, :Moderna=>61:66, :JnJ=>67:72, :totvaccinated=>73:78,            # vaccine
        :base=>79:84, :alpha=>85:90, :delta=>91:96, :omicron=>97:102                    # variant
        )   
    group = collect(keys(map2dict))  

    tmpdict = Dict{Int, Series}()   # dict of locales to Series

    for loc in locales
        tmpdict[loc] = Series(
            cum = zeros(Int, n_days, map2dict[last(group)][end]), 
            new = zeros(Int, n_days, map2dict[last(group)][end]), 
            groups = group,
            cols = map2dict
            )
    end

    return tmpdict       
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
    loadvariants = keys(infectdict) 

    for variant in loadvariants
        newdict = merge(infectdict[variant][:spread], infectdict[variant][:immunity])
        # newdict = Dict(Symbol(k) => v for (k,v) in newdict)
        infectset[Symbol(variant)] = Infectparams(newdict)
    end

    if isempty(variantlist)
        for variant in loadvariants
            push!(variantlist, Symbol(variant))   # this is a module global variable. Forgive me for I have sinned--except it makes sense...
        end
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
# Data mapping: for types and values
#####################################################################################

function mapcontact(x::condition)
    Int(x)-4
end

function mapage(x::agegrp)
    Int(x)
end

function maptouch(x::Union{condition, status})
    if x == unexposed
        1
    elseif x == recovered
        2
    elseif x == nil
        3
    elseif x == mild
        4
    elseif x == sick
        5
    elseif x == severe
        6
    else
        @assert false "invalid index to touchfactors $x"
    end
end

function tup2vec(maptup, vals)
    [getfield(maptup, x) for x in vals]
end


#= 
lookup tables for enum values: 
- don't need lookup for Int or Symbol: just use Symbol(nil) and Int(nil)-->these are faster than any lookup
- for symbol use symcond[:nil] => nil::condition = 5
- for string use symcond[Symbol("nil")] => nil::condition = 5
=#

"""
    symboltoagegrp(x::Union{Symbol, String})
Lookup a string or symbol that matches an enum value of Enum agegrp.
Generates an error if the string or symbol does not match.

Examples:
- symbol2agegrp(:age0_19) result:  agegrp::age0_19 = 1
- symbol2agegrp("age0_19") result: agegrp::age0_19 = 1
    
"""
function symbol2agegrp(x::Union{Symbol, String})::agegrp
    x = Symbol(x)
    inst_a = instances(agegrp)
    symtoage = freeze(Dict(zip(Symbol.(inst_a), inst_a))) # .5x time of regular dict
    @assert in(x, keys(symtoage)) "Error: input symbol $x is not an agegrp value."

    symtoage[x]
end


"""
    symbol2condition(x::Union{Symbol String})::condition  
Lookup a string or symbol that matches an enum value of Enum condition.
Generates an error if the string or symbol does not match.

Examples:
- symbol2cond(:nil)  result: nil::condition = 5
- symbol2cond("nil") result: nil::condition = 5
    
"""
function symbol2condition(x::Union{Symbol, String})::condition  
    x = Symbol(x) 
    inst_cond = instances(condition)
    symtocond = freeze(Dict(zip(Symbol.(inst_cond), inst_cond)))

    symtocond[x]
end


"""
    symbol2status(x::Union{Symbol String})::status  
Lookup a string or symbol that matches an enum value of Enum status.
Generates an error if the string or symbol does not match.

Examples:
- symbol2cond(:recovered)  result: recovered::status = 3
- symbol2cond("recovered") result: recovered::status = 3
    
"""
function symbol2status(x::Union{Symbol, String})::status  
    x = Symbol(x) 
    inst_status = instances(status)
    symtostatus = freeze(Dict(zip(Symbol.(inst_status), inst_status)))

    symtostatus[x]
end


"""
    symbol2allconds(x::Union{Symbol String})::Union{condition, status}  
Lookup a string or symbol that matches an enum value of Enum condition or status.
Generates an error if the string or symbol does not match.

Examples:
- symbol2allconds(:nil)  result: nil::condition = 5
- symbol2allconds("nil") result: nil::condition = 5
- symbol2allconds(:dead) result: dead::status = 4
    
"""
function symbol2allconds(x::Union{Symbol, String})::Union{condition, status}  
    x = Symbol(x) 
    inst_status = instances(status)
    inst_cond = instances(condition)

    symtostatus = freeze(Dict(zip(Symbol.(inst_status), inst_status)))
    symtocond = freeze(Dict(zip(Symbol.(inst_cond), inst_cond)))
    symtoallconds = merge(symtostatus, symtocond)

    symtoallconds[x]
end



# lookup table for shift
inst_shift = instances(shift)

"""
    symtoshift[sh::Symbol]
Dict used as lookup table to convert symbol or string to enum value for an shift.

Examples:
- for symbol use symtoshift[:recover] returns shift::recover = 1
- for string use symtoshift[Symbol("recover")] returns shift::recover = 1
    
"""
const symtoshift = freeze(Dict(zip(Symbol.(inst_shift), inst_shift))) # .5x time of regular dict


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


######################################################################################
# precalculate agegrp indices--these do not change during the simulation
######################################################################################


function precalc_agegrp_filt(dat)  # dat for a single locale
    agegrp_filt_bit = Dict(age => dat.agegrp .== age for age in agegrps)
    agegrp_filt_idx = Dict(age => findall(agegrp_filt_bit[age]) for age in agegrps)
    return (boolean=agegrp_filt_bit, idx=agegrp_filt_idx)
end
