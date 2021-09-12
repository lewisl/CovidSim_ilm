######################################################################################
# setup and initialization functions: ILM Model
######################################################################################


function setup(n_days, locales;  # must provide following inputs
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
        datadict = build_data(locales, geodata, n_days)

    # social parameters
        socialparams = build_socialparams(socialfilename, paramdir)

    # variants for spread parameters and transition decision trees
        variants = YAML.load_file(joinpath(paramdir, variantsfilename); dicttype=Dict{Symbol,Any})

    # spread parameters
        spreadset = build_spread_params(variants, paramdir)

    # transition arrays 
        (transitionset, trvec) = build_transition_params(variants, paramdir)

    # vaccines  TODO this is not the right approach
    if dovax
        vaxset = build_vaxset(vaccinefilename, paramdir=paramdir)
        vxschedset = build_vaxschedset()
    else
        vaxset = nothing
        vxschedset = nothing
    end

    return (dat=datadict, transitionset=transitionset, geo=geodata, vaxset=vaxset,
            vxschedset=vxschedset, spreadset=spreadset, social=socialparams, trvec = trvec)  
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
    # return Dict("popdat"=>popdat, "isolatedmx"=>isolatedmx, "testmx"=>testmx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
    return Dict("popdat"=>popdat, "agegrp_idx"=>agegrp_idx, "cumhistmx"=>cumhistmx, "newhistmx"=>newhistmx)
end


"""
Pre-allocate and initialize population data for one locale in the simulation.
"""
function pop_data(pop; age_dist=age_dist)

        parts = apportion(pop, age_dist)
        dat = Table(
            pid = collect(1:pop),  # ordinal persistent id for persons in matrix
            status = fill(unexposed, pop),    
            agegrp = reduce(vcat,[fill(age, parts[Int(age)]) for age in agegrps]), 
            cond = fill(notsick, pop),
            sickday = zeros(Int, pop),   
            variant = fill(:default, pop),
            recovday = zeros(Int, pop),  
            deadday = zeros(Int, pop),   
            cluster = zeros(Int, pop), 
            sdcomply = fill(:none, pop),  
            vaxstatus = fill(:none, pop),  # :none, :first, :multiple, :full, :booster  maybe others later...
            vaxrcvd = fill([:none], pop),    # vaccine symbols  :pfizer, :moderna, :jnj
            vaxday = fill([0], pop), 
            fullvaxday = zeros(Int, pop),
            test = falses(pop),  
            testday = zeros(Int, pop),  
            quar = falses(pop),
            quarday = zeros(Int, pop))

    return dat       
end


function hist_dict(locales, n_days; conds=allconds, agegrps=n_agegrps)
    dat = Dict{Int64, Array{Int}}()
    for loc in locales
        dat[loc] = zeros(Int, n_days, last(last(map2series))) # (conds, agegrps + 1, n_days) => (8, 6, 150)
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
    spreadset = Dict{Symbol, Union{Infectparams, Vaccineparams}}()
    for variant in keys(variants)
        v = YAML.load_file(joinpath(paramdir, "variant_parameters", variants[variant][:directory_name],
            variants[variant][:infect_fname]), dicttype=Dict{Symbol, Any})
        v = Infectparams(v)
        spreadset[variant] = v
    end
    return spreadset
end


function build_transition_params(variants, paramdir)
    transitionset = Dict()

    for variant in keys(variants)
        transitionset[variant] = setup_dt(joinpath(paramdir, "variant_parameters", variants[variant][:directory_name], 
            variants[variant][:transition_fname])) 
    end

    # pre-allocate trvec used in hot loop: no. of columns in transition array
    sz = size(first(first(transitionset[:default])[end])[end][:transition], 2)
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
