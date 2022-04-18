
# an approach to callback functions
# NONE OF THIS WORKS.  ALL RESULT IN CLOSURES

function passtheargs(;a, b, c)
       (x=a, y=b, z=c)
end

function addup(; x,y,z)
    x + y + z
end

# val10=5; val20=7; val30=11  # block closure by commenting this out


# another way

function blowup(;x=a, y=b, z=c)
    x*y.+z
end

a=3; b=1.5; c=10

blowup()  # should be 15

# avoid globals

function testit(a,b,c)
    val10 = a + 2
    val20 = b * 2
    val30 = c .* 0.5
    func = blowup
    argset = passtheargs(a=val10, b=val20, c=val30)
    @show Base.summarysize(argset)
    b1 = func(;argset...)
    
end