


@Base.kwdef mutable struct Vaccineparams
    name::Symbol = :Pfizer
    shots::Int
    halflife::Int
    sendrisk::Vector{Float64}
    recvrisk::Vector{Float64}
    recvrisk_reduction::Dict{Symbol, Float64}
    pattern::Vector{Float64}
end

"""
Method for converting a dict created from YAML to this struct
"""
Vaccineparams(vd) =
    (Vaccineparams(
        name                     = Symbol(vd[:name]),
        shots                    = vd[:shots],
        halflife                 = vd[:halflife],
        sendrisk                 = vd[:sendrisk],
        recvrisk                 = vd[:recvrisk],
        recvrisk_reduction       = vd[:recvrisk_reduction], 
        pattern                  = vd[:pattern]
    ))


@Base.kwdef mutable struct Vaxsched
    vaccine::Symbol # use this as the name of the schedule, not the name of the vaccine.
    dayrange::UnitRange{Int64}
    targetpct::Float64
    pctperdayfn::Function
end


"""
    makevaxfn(dayrange, pattern)

    returns: function pctperday(day)

Create a function that implements the vaccination schedule for a given vaccine type.
This function, when called with a simulation day, returns the percentage of the
target population to receive the vaccine on that day.
"""
function makevaxfn(dayrange, pattern, targetpct)
    schedlength = length(dayrange)
    distribscale = pattern .* ((length(pattern)-1)/schedlength)
    interp = LinearInterpolation(0:length(distribscale)-1, distribscale)
    startday = dayrange.start; endday = dayrange.stop
    return  function pctperday(day)
                @assert startday <= day <= endday "Day not in dayrange $dayrange"
                p = (day - startday + 1) * round((length(interp)-1) / schedlength, digits=4)
                p = p < length(pattern) - 1 ? p : length(pattern) - 1 
                return interp(p) * targetpct
            end
end



function setupvax()
end


function getashot!(locdat, p, day, vx, vaxkeys)
    @assert (vx in vaxkeys) "$vx is not a valid vaccine in $vaxkeys"
    if isnothing(locdat.vax[p])
        locdat.vax[p] = [vx]
        locdat.vaxday[p] = [day]
    else
        push!(locdat.vax[p], vx)
        push!(locdat.vaxday[p], day)
    end
end


function getashot!(locdat, pvec::Union{Vector{Int}, UnitRange{Int}}, day, vx, vaxkeys)
    for p in pvec
        getashot!(locdat, p, day, vx, vaxkeys)
    end
end





