function getval(outer, s1, s2)
    getproperty(getproperty(outer, s1), s2)
end

function getkey(outer, s1, s2)
    get(get(outer, s1, nothing), s2, nothing)
end


nt = @benchmarkable  getval(b, g2, g3) setup = begin; 
                                    b = (foo=(a=1, b=2), bar=(c=3, d=4)); 
                                    g2 = :foo; g3 = :a; 
                                end;

dct = @benchmarkable getkey(b, g2, g3) setup = begin; 
                                    b = Dict(:foo => Dict(:a=>1, :b=>2), :bar => Dict(:c=>3, :d=>4)); 
                                    g2=:foo; 
                                    g3=:bar; 
                                end;