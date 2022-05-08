# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %%
n1 = 100_000
#poparr = Array{Union{Int8, Int16},2}(undef, n1, 6)

poparr = zeros(Int16,n1,6)
println(Base.summarysize(poparr))

# %%
age_dist = [0.251, 0.271,   0.255,   0.184,   0.039]

# %%
using Distributions
using IndexedTables
using CategoricalArrays

# %%
dcat = Categorical(age_dist)

# %%
poparr[:,1] = rand(dcat,n1);
poparr[:,3] = ones(n1);
@time sort!(poparr,dims=1);

# %%
@time [count.(poparr[:,1] .== i for i in 1:5)] ./ n1

# %%
println(typeof(poparr))
Base.summarysize(poparr)

# %%
n = 500_000
poptable = table(zeros(Int8,n),zeros(Int8,n),zeros(Int8,n),zeros(Int16,n),zeros(Int16,n),zeros(Int8,n),
    names=[:agegrp,:sick,:status,:tested,:isolated,:duration], pkey=[:agegrp,:sick,:status,:tested,:isolated,:duration]);

# %%
Base.summarysize(poptable)

# %%
select(poptable, 1) .= Int8.(rand(dcat,n));
select(poptable, 3) .= ones(Int8, n);

# %%
@time poptable = table(sort(poptable, :agegrp), pkey=[:agegrp,:sick,:status,:tested,:isolated,:duration])

# %%
@time select(groupby(length, poptable, :agegrp), :length) ./ n

# %%
filter(r->r.agegrp == 5, poptable);

# %%
