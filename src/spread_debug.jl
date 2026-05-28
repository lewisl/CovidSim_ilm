#########################################################################################
# spread_debug.jl: event-level spread tracing for Julia/C++ parity checks
#########################################################################################

Base.@kwdef struct SpreadDebugConfig
    max_days::Int = 3
    max_spreaders::Int = 250
    max_contacts::Int = 2000
end

function init_spread_debug_trace(locales; config=SpreadDebugConfig())
    Dict(loc => (
            config = config,
            spreaders = (
                day = Int[],
                spreader_id = Int[],
                spr_agegrp = Symbol[],
                spr_cond = Symbol[],
                spr_duration = Int[],
                spr_variant = Symbol[],
                indoor_factor = Float64[],
                density_factor = Float64[],
                contact_factor = Float64[],
                contact_scale = Float64[],
                num_contacts = Int[],
                sendrisk = Float64[],
            ),
            contacts = (
                day = Int[],
                spreader_id = Int[],
                contact_order = Int[],
                contact_id = Int[],
                targ_agegrp = Symbol[],
                targ_status = Symbol[],
                targ_cond = Symbol[],
                indoor_factor = Float64[],
                touch_factor = Float64[],
                touch_prob = Float64[],
                touched = Bool[],
                sendrisk = Float64[],
                recvrisk = Float64[],
                recovfactor = Float64[],
                vaxfactor = Float64[],
                infect_risk = Float64[],
                infected = Bool[],
            ),
        ) for loc in locales)
end

@inline spread_debug_enabled(trace) = !isnothing(trace)
@inline spread_debug_rowcount(cols) = length(cols.day)

@inline function should_trace_spreader(trace, thisday)
    spread_debug_enabled(trace) || return false
    thisday <= trace.config.max_days || return false
    spread_debug_rowcount(trace.spreaders) < trace.config.max_spreaders
end

@inline function should_trace_contact(trace, thisday)
    spread_debug_enabled(trace) || return false
    thisday <= trace.config.max_days || return false
    spread_debug_rowcount(trace.contacts) < trace.config.max_contacts
end

function push_spreader_debug!(trace; day, spreader_id, spr_agegrp, spr_cond, spr_duration, spr_variant,
        indoor_factor, density_factor, contact_factor, contact_scale, num_contacts, sendrisk)
    cols = trace.spreaders
    push!(cols.day, day)
    push!(cols.spreader_id, spreader_id)
    push!(cols.spr_agegrp, spr_agegrp)
    push!(cols.spr_cond, spr_cond)
    push!(cols.spr_duration, spr_duration)
    push!(cols.spr_variant, spr_variant)
    push!(cols.indoor_factor, indoor_factor)
    push!(cols.density_factor, density_factor)
    push!(cols.contact_factor, contact_factor)
    push!(cols.contact_scale, contact_scale)
    push!(cols.num_contacts, num_contacts)
    push!(cols.sendrisk, sendrisk)
    return trace
end

function push_contact_debug!(trace; day, spreader_id, contact_order, contact_id, targ_agegrp, targ_status,
        targ_cond, indoor_factor, touch_factor, touch_prob, touched, sendrisk, recvrisk, recovfactor,
        vaxfactor, infect_risk, infected)
    cols = trace.contacts
    push!(cols.day, day)
    push!(cols.spreader_id, spreader_id)
    push!(cols.contact_order, contact_order)
    push!(cols.contact_id, contact_id)
    push!(cols.targ_agegrp, targ_agegrp)
    push!(cols.targ_status, targ_status)
    push!(cols.targ_cond, targ_cond)
    push!(cols.indoor_factor, indoor_factor)
    push!(cols.touch_factor, touch_factor)
    push!(cols.touch_prob, touch_prob)
    push!(cols.touched, touched)
    push!(cols.sendrisk, sendrisk)
    push!(cols.recvrisk, recvrisk)
    push!(cols.recovfactor, recovfactor)
    push!(cols.vaxfactor, vaxfactor)
    push!(cols.infect_risk, infect_risk)
    push!(cols.infected, infected)
    return trace
end

@inline spread_debug_tables(trace) = (
    spreaders = LazyTable(; (name => value for (name, value) in pairs(trace.spreaders))...),
    contacts = LazyTable(; (name => value for (name, value) in pairs(trace.contacts))...),
)
