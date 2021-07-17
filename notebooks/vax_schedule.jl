# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %%
using Interpolations
using Polynomials
using Plots

# %%
@Base.kwdef struct Vaxsched
    dayrange::UnitRange{Int64}
    targetpct::Float64
    shape::Vector{Float64}
end

# %%
function makevax(vx::Vaxsched)
    schedlength = length(vx.dayrange)
    shapescale = vx.shape .* ((length(vx.shape)-1)/schedlength)
    interp = LinearInterpolation(0:length(shapescale)-1, shapescale)
    startday = vx.dayrange.start; endday = vx.dayrange.stop
    return  function vaxsched(day)
                @assert startday <= day <= endday "Day not in dayrange $(vx.dayrange)"
                p = (day - startday + 1) * round((length(interp)-1) / schedlength, digits=4)
                return interp(p) * vx.targetpct
            end
end

# %% tags=[]
maxpct = .85
shape = [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01]
@show sum(shape), length(shape)
xs = 91:360
id = collect(0:10)
order = 10
shapescale = shape .* ((length(shape)-1)/length(xs))

# %% tags=[]
vx = Vaxsched(dayrange=xs, targetpct=maxpct, shape=shape)
dump(vx)

# %% [markdown]
# ## Interpolate with polynomial fit

# %%
px = fit(collect(0:length(shapescale)-1), shapescale, order)

# %%
function relx(p, shape, xs) 
    startday = first(xs); endday = last(xs)
    @assert startday <= p <= endday "Day not in startday to endday"
    (p - startday + 1) * round(((length(shape)-1) / length(xs)), digits=4)
end

# %%
relx(101, shapescale, xs)

# %%
11/270

# %%
allx = [relx(i, shape, xs) for i in xs]

# %%
py = px.(allx)

# %%
py .= max.(0.0, py)

# %%
sum(py)

# %%
# This is too artificial
plot(allx, py)

# %%
@show allx[length(xs)] py[length(xs)]

# %%
sum(py .* .65 .* 330_000_000)

# %%
327_000_000

# %% [markdown]
# ## piecewise linear

# %%
interp = LinearInterpolation(0:length(shapescale)-1, shapescale)

# %%
interp.(allx)[1]

# %%
vaxsched = makevax(vx)

# %%
vaxsched(91)

# %%
plot(interp(allx), size=(800,400))

# %%
plot(xs,vaxsched.(xs), size=(800,400))

# %%
sum(vaxsched.(xs))

# %%
