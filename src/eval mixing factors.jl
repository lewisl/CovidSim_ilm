using TypedTables
using Plots
alts = Table(fac1=0.03:0.01:1.0, fac2=0.03:0.01:1.0)
alts = @Select(fac1, fac2, fac3=reverse($fac2), prod1=$fac1 .* $fac2)(alts)
alts = @Select(fac1, fac2, fac3, prod1, prod2=$fac1 .* $fac3, min=min.($fac1, $fac3))(alts)
alts = @Select(fac1, fac2, fac3, prod1, prod2, min1=$min, exp1=$prod2 .* exp.(0.85 .- $prod2))(alts)

alts = Table(alts, @Select(exp2=$prod2 .* exp.(0.2 .- $prod2))(alts))

plot(alts.min1)
plot!(alts.exp1)
plot!(alts.exp2)

# Squashing functions
xs = range(0.01,3.0,300)

