module ThermoCycleBench
    using BenchmarkTools, Clapeyron, ThermoCycleGlides
    
    const p = 101325.0
    const h = collect(range(-1,1,100));
    const fluid  = cPR(["Propane","R134a"],idealmodel = ReidIdeal); 
    const z = [1.0,1.0]
    const p_out = 5*p
    const η_isen = 0.8
    f_compression(x) = ThermoCycleGlides.isentropic_compressor(p,p_out,η_isen,x,z,fluid)
    suite = BenchmarkGroup()
    suite["utils"] = BenchmarkGroup()
    suite["utils"]["compressor"] = @benchmarkable f_compression.(x) setup=(x = h)
end

ThermoCycleBench.suite

