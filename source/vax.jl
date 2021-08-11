
###################################################
# data structures for vaccination
###################################################

@Base.kwdef mutable struct Vaccineparams
    name::Symbol 
    shots::Int
    delay2ndshot::Union{Int, Nothing}   # days until 2nd shot (probability less important)
    halflife::Int  # days to 50% decline in effectiveness
    sendrisk::Vector{Float64}
    recvrisk::Vector{Float64}
    recvrisk_reduction::Dict{Symbol, Float64}
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        Vaccineparams(vd::Dict)=(
            Vaccineparams(
                name                     = Symbol(vd[:name]),
                shots                    = vd[:shots],
                delay2ndshot             = vd[:delay2ndshot],
                halflife                 = vd[:halflife],
                sendrisk                 = vd[:sendrisk],
                recvrisk                 = vd[:recvrisk],
                recvrisk_reduction       = vd[:recvrisk_reduction], 
                )
        )



@Base.kwdef mutable struct Vaxsched
    vax::Dict{Symbol, NamedTuple{(:mix, :doses, :pct2ndshot), Tuple{Float64, Int64, Float64}}}
    dayrange::UnitRange{Int64}
    targetpct::Float64
    filterfunc::Function  # discretionary criteria for who is included in the vaccine schedule
    shotmode::Symbol   # one of :first, :second, :all.  best to use :all. others force only single shot
    # pattern defines 11 point on a piecewise linear "curve" of the pace of vaccination as the
    # percentage of a given population that receives vaccination-> must sum to 1.0
    # [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]
    pattern::Vector{Float64} # = [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]         
    pctfunc::Function
end


# this works
"""
Method to build schedule of people to be vaccinated by percent by day.
This will typically be done as part of a runcase rather than
as part of initial setup. ???
"""
function create_vaxsched(;vax::Dict{Symbol, NamedTuple{(:mix, :doses, :pct2ndshot), Tuple{Float64, Int64, Float64}}},
                        dayrange, targetpct, 
                        filterfunc=basevaxfilterfunc,
                        pattern=[0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01], 
                        shotmode=:all)

    vaxmix = [v.mix for v in values(vax)]
    @assert sum(vaxmix) == 1.0 "Sum of mixes for vaccines does not equal 1: $vaxmix"

    Vaxsched(
        vax=vax,
        dayrange=dayrange,
        targetpct=targetpct,
        filterfunc=filterfunc,
        shotmode=shotmode,
        pattern=pattern,
        pctfunc=makevaxpctperdayfn(dayrange, targetpct, pattern; shotmode=:all)
        )
end


###########################################################
# setup code for vaccination
###########################################################

# this works
function setupvax(paramdir="../parameters")
    vaccinefiles = YAML.load_file(joinpath(paramdir, "vaccines.yml"); dicttype=Dict{Symbol,Any})

    pprint(vaccinefiles)

    vaxset = build_vaxset(vaccinefiles, paramdir)

    return vaxset
end


# not used yet
function setupvaxsched(vaccines, paramdir="../parameters")
    vxschedset = build_vxschedset(vaccines, paramdir)
end


"""
    makevaxfn(dayrange, pattern)

    returns: function pctperday(day)

Create a function that implements the vaccination schedule for a given vaccine type.
This function, when called with a simulation day, returns the percentage of the
target population to receive the vaccine on that day.
"""
function makevaxpctperdayfn(dayrange, targetpct, pattern; shotmode=:all)
    schedlength = length(dayrange)
    distribscale = pattern .* ((length(pattern)-1)/schedlength)
    interp = LinearInterpolation(0:length(distribscale)-1, distribscale)
    startday = dayrange.start; endday1 = dayrange.stop
    function vaxpctperday(day)
        p = (day - startday + 1) * round((length(interp)-1) / schedlength, digits=4)
        p = p < length(pattern) - 1 ? p : length(pattern) - 1 
        return interp(p) * targetpct
    end
end

# this works
function build_vaxset(vaccinefiles, paramdir="../parameters")
    vaxset = Dict{Symbol, Dict{Symbol, Union{Vaccineparams, transitionT}}}()
    for vax in keys(vaccinefiles)
        vparams = YAML.load_file(joinpath(paramdir,vaccinefiles[vax][:directory_name],
            vaccinefiles[vax][:infect_fname]), dicttype=Dict{Symbol, Any})
        vparams = Vaccineparams(vparams)
        vtrans = setup_dt(joinpath("../parameters", vaccinefiles[vax][:directory_name], 
                    vaccinefiles[vax][:transition_fname]))
        vaxset[vax] = Dict()
        vaxset[vax][:params] = vparams
        vaxset[vax][:transition] = vtrans
    end
    return vaxset
end


# function build_vax_transition(vaccines, paramdir="../parameters")
#     vaxtranset = Dict()
#     for vax in keys(vaccines)
#         vaxtranset[vax] = setup_dt(joinpath("../parameters", vaccines[vax][:directory_name], 
#                 vaccines[vax][:transition_fname]))
#     end
#     return vaxtranset
# end



# method for one person
function giveshot!(locdat, p, day, vx, shotmode)
    if shotmode == :1st
        locdat.vax[p] = [vx]
        locdat.vaxday[p] = [day]
    else
        push!(locdat.vax[p], vx)
        push!(locdat.vaxday[p], day)
    end
end

# method for multiple people--not sure we'll ever use given loop in vaccinate!
function giveshot!(locdat, pvec::Union{Vector{Int}, UnitRange{Int}}, day, vx, shotmode)
    for p in pvec
        getashot!(locdat, p, day, vx, shotmode)
    end
end


function vaccinate!(locdat, vxschedset, contactable_idx, spreadset, dovax, dovariant)
    
    dovax || return
    
    today = day_ctr[:day]

    for vxsched in vxschedset

        # shortcircuit the whole shebang
        dayrange = vxsched.dayrange
        delay2ndshot = Dict(v=>spreadset[v].delay2ndshot for v in vaccines)
        maxstop2ndshot = dayrange.stop + maximum(values(delay2ndshot))

        if (today < dayrange.start) | (today > maxstop2ndshot)
            continue
        end

        # schedule parameters
        filterfunc = vxsched.filterfunc
        shotmode = vxsched.shotmode # values in :first, :second, :all
        pctfunc = vxsched.pctfunc
        vaccines = keys(vxsched.vaccines)
        vaxprops = vxsched.vax

        # vaccine parameters:  vaccine=>trait
        shots = Dict(v=>spreadset[v].shots for v in vaccines)
        # see delay2ndshot above...

        # day flags
        giveshot1 = startvax <= today <= stop1stshot

        # counters
        shot1todaycnt = 0
        shot2todaycnt = 0
        # NOT RIGHT ALSO Doses maxshot1s = vaxlpctperdayfn(day) * length(tovaxidx)

        # columns
        vaxstatuscol = locdat.vaxstatus
        vaxdaycol = locdat.vaxday
        shotscol = locdat.shots
        statuscol = locdat.status
        condcol = locdat.cond
        agegrpcol = locdat.agegrp


        if giveshot1

        else 
            # second+ shots

        end
    end  # for vxsched

end


function basevaxfilterfunc(p, agegrpcol)
    (agegrpcol[p] != age0_19)
end