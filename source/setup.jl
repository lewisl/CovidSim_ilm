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
    variantsfilename)

    # geodata
        geodata = buildgeodata(geofilename)

    # simulation data matrix
        datadict = build_data(locales, geodata, ndays)

    # social parameters
        socialparams = build_socialparams(socialfilename, paramdir)

    # variants for spread parameters and transition decision trees
        variants = YAML.load_file(joinpath(paramdir, variantsfilename); dicttype=Dict{Symbol,Any})

    # spread parameters
        spreadset = build_spread_params(variants, paramdir)

    # transition arrays 
        (transitionset, trvec) = build_transition_params(variants, paramdir)

    # vaccines  TODO this is not the right approach: test if we have vax inputs instead
    if dovax
        vaxset = build_vaxset(vaccinefilename, paramdir=paramdir)
        vaxschedset = build_vaxschedset()
    else
        vaxset = nothing
        vaxschedset = nothing
    end

    return (ndays = ndays, locales=locales, dat=datadict, transitionset=transitionset, geo=geodata, vaxset=vaxset, variants=variants,
            vaxschedset=vaxschedset, spreadset=spreadset, social=socialparams, trvec = trvec)  
end


"""
Convert a vector of dates from a csv file in format "mm/dd/yyyy"
to a vector of Julia numeric Date values in format yyyy-mm-dd
"""
function quickdate(strdates)  # 20x faster than the built-in date parsing, e.g.--runs in 5% the time
    ret = [parse.(Int,i) for i in split.(strdates, '/')]
    ret = [Date.(i[3], i[1], i[2]) for i in ret]
end


function build_data(locales, geodata, n_days)

    pop = [geodata[geodata[:, "fips"] .== loc, "pop"][1] for loc in locales]

    popdat = Dict(loc => pop_data(geodata[geodata[:, "fips"] .== loc, "pop"][1]) for loc in locales)

    # precalculate agegrp indices
    agegrp_idx = Dict(loc => precalc_agegrp_filt(popdat[loc]).idx for loc in locales)

    cumhistmx = hist_dict(locales, n_days)
    newhistmx = hist_dict(locales, n_days)
    
    return Dict("popdat"=>popdat, "agegrp_idx"=>agegrp_idx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
end


"""
Pre-allocate and initialize population data for one locale in the simulation.
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
const map2series = (unexposed=1:6, infectious=7:12, recovered=13:18, dead=19:24,          # status
                    nil=25:30, mild=31:36, sick=37:42, severe=43:48, totinfected=49:54,   # conditions
                    Pfizer=55:60, Moderna=61:66, JnJ=67:72, totvaccinated=73:78,          # vaccines
                    base=79:84, alpha=85:90, delta=91:96, omicron=97:102)                 # variants


function hist_dict(locales, n_days; conds=allconds, agegrps=n_agegrps)
    dat = Dict{Int64, Array{Int}}()
    for loc in locales
        dat[loc] = zeros(Int, n_days, map2series[end][end]) 
    end
    return dat       
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


function build_spread_params(variants, paramdir)
    spreadset = Dict{Symbol, Infectparams}()
    for variant in keys(variants)
        v = YAML.load_file(joinpath(paramdir, "variant_parameters", variants[variant][:directory_name],
            variants[variant][:infect_fname]), dicttype=Dict{Symbol, Any})
        v = Infectparams(v)
        spreadset[variant] = v  # access a param as spreadset[:alpha].recvrisk
    end
    return spreadset
end


function build_infect_params(variantfilename, paramdir)
    infectdict = YAML.load_file(joinpath(paramdir, variantfilename), dicttype=Dict{Symbol, Any})

    
end


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


function build_socialparams(socialfilename, paramdir)

    social_inputs = YAML.load_file(joinpath(paramdir, socialfilename))

    required_params = ["contactfactors", "touchfactors", "gammashape"]
    has_all = true
    lacking = []
    for p in required_params
        if !haskey(social_inputs, p)
            push!(lacking, p)
            has_all = false
        end
    end
    @assert has_all "required keys: $lacking not in $(infectfilename)"

    Socialparams(
        gammashape         = social_inputs["gammashape"],
        contactfactors    = Dict(symtoage[Symbol(k1)] => 
                                Dict(symtocond[Symbol(k2)] => Float64(v2) for (k2, v2) in v1)  for (k1, v1) in social_inputs["contactfactors"]),
        touchfactors      = Dict(symtoage[Symbol(k1)] => 
                                Dict(symtoallconds[Symbol(k2)] => Float64(v2) for (k2, v2) in v1)  for (k1, v1) in social_inputs["touchfactors"])
        )
    
end



#####################################################################################
# helper functions for setup
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


function map2vec(maptup, vals)
    [getfield(maptup, x) for x in vals]
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
# agegrp_filt_bit, agegrp_filt_idx = precalc_agegrp_filt(ilmat);
