abstract type Vaccine end

    @Base.kwdef struct Pfizer <: Vaccine
        name::Symbol
        shots::Int
        halflife::Int
        send_risk::Vector{Float64}
        recv_risk::Vector{Float64}
    end

@Base.kwdef struct Vaxsched
    dayrange::UnitRange{Int64}
    targetpct::Float64
    shape::Vector{Float64}
end


function makevax(vx::Vaxsched)
    schedlength = length(vx.dayrange)
    shapescale = vx.shape .* ((length(vx.shape)-1)/schedlength)
    interp = LinearInterpolation(0:length(shapescale)-1, shapescale)
    startday = vx.dayrange.start; endday = vx.dayrange.stop
    return  function vaxsched(day)
                @assert startday <= day <= endday "Day not in dayrange $(vx.dayrange)"
                p = (day - startday + 1) * round((length(interp)-1) / schedlength, digits=4)
                return interp(p) * maxpct
            end
end



function setup_vax()
end

function dovax()
end

    """
    isinfected(riskmx, spreadersickday, contactagegrp)::Bool

Returns true if the spreader infected the contact. 
"""
@inline function isinfected(vax_spreadparams, spreadparams,  spreader, contact, locdat)::Bool
    @inbounds @fastmath prob = (spreadparams.send_risk[spreadersickday] * 
                        spreadparams.recv_risk[Int(contactagegrp)])            # TODO also vaccinated people will have partially unsusceptible
    return @fastmath rand(Binomial(1, prob)) == 1
end

function make_isinfected(spreadparams)
    return 
        
end


function shots(vaxname::Symbol, vaxsched)