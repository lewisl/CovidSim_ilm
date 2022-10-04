# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: julia 1.8 threads auto env v1.8 1.8.0
#     language: julia
#     name: julia-1.8-threads-auto-env-v1.8-1.8
# ---

# %%
using DataFrames
using Distributions
using Plots
using Random
using StatsBase
using Markdown

# %% [markdown]
# ## Figuring out what distribution and paramaters to use; how to apply to events
#
# ### Binomial

# %%
dbin = Binomial(1000,0.1)
x = rand(dbin,100);
println(mean(x))
println(x)
histogram(x, bins=10)


# %%
# bin the data
histdata = [count(isequal(i),x) for i in unique(x)]
println(histdata)
# hh = hist(x,22);  # using PyPlot package
# println(hh[1])
# println(hh[2])


# %%
dhypergeo = Hypergeometric(100,200,100) # successes, failures, draws
x = rand(dhypergeo,100)
println(mean(x))
histogram(x,bins=10)

# %%
dpoiss = Poisson(1)
x = rand(dpoiss, 1000)
println(mean(x))
histogram(x)

# %% [markdown]
# ### The long tail of (some) gamma distributions.

# %%
dgamma = Gamma(1.4, 2.0 )  #shape, scale
x = rand(dgamma,100);
println(mean(x))
println(round(sum((x))))
histogram(x, bins=22)
# hh = hist(x,20);  using PyPlot
# println(ceil(sum(hh[1] .* hh[2][2:end])))
# println(mean(dgamma) * 20)

# %%
derlang = Erlang(1,5)  # integer shape, scale
x = rand(derlang,100);
histogram(x,bins=22)


# %%
println(round.(Int,x))
println(round(mean(x)))


# %% [markdown]
# ### Variations in Gamma Parameters

# %%
# multiple gamma distributions, varying parameters shape, scale
shapes = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]
scales = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5]
p = Array{Any,2}(undef,6,6)
for (shi,sh) in enumerate(shapes)
    for (sci, sc) in enumerate(scales)
        Random.seed!(1)
        dgamma = Gamma(sh, sc)
        x = rand(dgamma,100)
        p[sci,shi] = histogram(x,bins=11,title="Shape $sh, Scale $sc", 
                               titlefontsize=6, size=(800,500))
    end
end
l = @layout([a1 a2 a3 a4 a5 a6; 
             b1 b2 b3 b4 b5 b6; 
             c1 c2 c3 c4 c5 c6; 
             d1 d2 d3 d4 d5 d6; 
             e1 e2 e3 e4 e5 e6; 
             f1 f2 f3 f4 f5 f6])
display(md"# Gamma Distribution Parameters")
display(md"Rows: Shape        Columns: Scale")
plot(p..., layout=l, legend=false, yaxis=false, tickfontsize=6)


# %% [markdown]
# ### Variations in Erlang Parameters

# %%
# multiple gamma distributions, varying parameters shape, scale
shapes = [1, 2, 3, 4, 5, 6]
scales = [1, 2, 3, 4, 5, 6]
p = Array{Any,2}(undef,6,6)
for (shi,sh) in enumerate(shapes)
    for (sci, sc) in enumerate(scales)
        Random.seed!(1)
        derlang = Erlang(sh, sc)
        x = rand(derlang,100)
        p[sci,shi] = histogram(x,bins=11,title="Shape $sh, Scale $sc", titlefontsize=6, size=(800,500))
    end
end
l = @layout([a1 a2 a3 a4 a5 a6; 
             b1 b2 b3 b4 b5 b6; 
             c1 c2 c3 c4 c5 c6; 
             d1 d2 d3 d4 d5 d6; 
             e1 e2 e3 e4 e5 e6; 
             f1 f2 f3 f4 f5 f6])

display(md"# Erlang Distribution Parameters")
plot(p..., layout=l, legend=false, yaxis=false, tickfontsize=6)


# %% [markdown]
# #### The pdf is a static characteristic of the distribution
# So it is always the same result.

# %%
dgamma = Gamma(1.2,1.6)
startpoint = 0.5
endpoint = 10
display(pdf.(dgamma,startpoint:endpoint))
println(sum(pdf.(dgamma,startpoint:endpoint)))

# %%
dgamma = Gamma(0.5,1.0)
x = rand(0.1:0.1:1.5,6)
rand(dgamma,6)

# %% jupyter={"outputs_hidden": true} tags=[]
methods(fit)

# %%
nilprobs = [.6, 0.0, .1, .1, 0.1, 0.1]

dcat = Categorical(nilprobs)
x = rand(dcat,100)
println(x)
println([count(x->x==i, x) for i in 1:6])
hh = fit(Histogram, x, nbins=6)
println(fieldnames(typeof(hh)))
println(hh.weights)
println(length(hh.weights))
println(hh.edges)

# %% [markdown]
# ## Categorical Distribution

# %% [markdown]
# ### Here is the categorical distribution output with 3 classes (well, 6 with 3 having 0.0 probability) with 100 trials. Next, we will try with similar inputs using the multinomial probability distribution.

# %%
snorm(arr) = arr ./ (sum(arr))

age_amplify = Dict(1=>(1.3,.7), 2=>(1.1,.9), 3=>(1.0, 1.0), 4=>(0.9,1.1), 5=>(.8,1.2))
nilprobs = [.6, 0.0, .3, .1, 0.0, 0.0]
# do we amplify the probs or the outcomes?
# probs
nilprobs = snorm(nilprobs .* [0.8, 0.0, 1.2, 1.2, 0.0, 0.0])
# all this does is jiggle things in a way that doesn't have a lot of meaning--it changes the relative probs in 
# partially obvious way.  We need a reason that justifies the math.

probs = nilprobs
println(probs)

dcat = Categorical(probs)
x = rand(dcat,100)
vals = [count(x .== i) for i in 1:6]
# names = ("cat 1", "cat 2", "cat 3", "cat 4", "cat 5", "cat 6")
# bar(names, vals);
println(x)
println(vals)
histogram(x,bins=10, xaxis=(lims=(1,10)))


# %%
# use the modified probs from the above to compare distribution plots
probs = [0.5, 0.0, 0.375, 0.125, 0.0, 0.0]
dmulti = Multinomial(100, probs)
names = ("cat 1", "cat 2", "cat 3", "cat 4", "cat 5", "cat 6")
x = vec(rand(dmulti, 1))
println(x)
bar(names, x);
# xlim(0.0,6.0)

# %% [markdown]
# ### Let's try this again with a simpler approach:

# %%
probs = [0.5, 0.0, 0.375, 0.125, 0.0, 0.0]
cats = 6
boundaries = [0.0, probs[1], sum(probs[1:2]), sum(probs[1:3]), sum(probs[1:4]), sum(probs[1:5]), sum(probs[1:6])]
println(boundaries)
x = rand(100)
counts = [count(boundaries[i-1] .<= x .< boundaries[i]) for i in 2:7]
names = ("cat 1", "cat 2", "cat 3", "cat 4", "cat 5", "cat 6")
println(counts)
bar(names, counts);

# %%
dbin1 = Binomial(100,0.5)
dbin3 = Binomial(100, 0.15)
dbin4 = Binomial(100, 0.05)

x1 = rand(dbin1, 100)
x3 = rand(dbin3, 100)
x4 = rand(dbin4, 100)

println(x1)
println(x3)
println(x4)

p1 = histogram(x1)
p3 = histogram(x3)
p4 = histogram(x4)

plot(p1,p3,p4,layout=@layout([a;b;c]),size=(400,700))

# %%
hist(rand(dbin1,1000));

# %%
dbin5 = Binomial(1,0.6)
x = rand(dbin5, 100)
println(sum(x))

# %%
