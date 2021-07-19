abstract type Vaccine end

    @Base.kwdef mutable struct Pfizer <: Vaccine
        name::Symbol = :Pfizer
        shots::Int
        halflife::Int
        send_risk::Vector{Float64}
        recv_risk::Vector{Float64}
    end

abstract type Vaxsched end

    @Base.kwdef mutable struct Pfizer_sched <: Vaxsched
        name::Symbol = :Pfizer # use this as the name of the schedule, not the name of the vaccine.
        dayrange::UnitRange{Int64}
        targetpct::Float64
        distrib::Vector{Float64} = [0.0, .02, .05, .10, .15, .19, 
                                    .21, .16, .08, .03, .01]
    end


"""
    makevax(vx::Vaxsched)

    returns: function vaxsched(day)

Create a function that implements the vaccination schedule for a given vaccine type.
This function, when called with a simulation day, returns the percentage of the
target population to receive the vaccine on that day.
"""
function makevax(vx::Vaxsched)
    schedlength = length(vx.dayrange)
    distribscale = vx.distrib .* ((length(vx.distrib)-1)/schedlength)
    interp = LinearInterpolation(0:length(distribscale)-1, distribscale)
    startday = vx.dayrange.start; endday = vx.dayrange.stop
    return  function vaxpctday(day)
                @assert startday <= day <= endday "Day not in dayrange $(vx.dayrange)"
                p = (day - startday + 1) * round((length(interp)-1) / schedlength, digits=4)
                return interp(p) * maxpct
            end
end



function setup_vax()
end


function getashot!(locdat, p, day, vx::Vaccine)
    if isnothing(locdat.vax[p])
        locdat.vax[p] = [vx.name]
        locdat.vaxday[p] = [day]
    else
        push!(locdat.vax[p], vx.name)
        push!(locdat.vaxday[p], day)
    end
end


function getashot!(locdat, pvec::Vector{Int}, day, vx::Vaccine)
    for p in pvec
        if isnothing(locdat.vax[p])
            locdat.vax[p] = [vx.name]
            locdat.vaxday[p] = [day]
        else
            push!(locdat.vax[p], vx.name)
            push!(locdat.vaxday[p], day)
        end
    end
end





