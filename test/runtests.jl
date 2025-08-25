using ThermoCycleGlides, Clapeyron
using Test


include("fluids.jl")

# @testset "ThermoCycleGlides.jl" begin
#     # Write your tests here.
# end


@testset "Isentopic Compression - Single Component" begin

    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]
        T1 = 300; p1 = 101325; p2 = 101325*5;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= ThermoCycleGlides.isentropic_compressor(p1,p2,1.0,h1,z1,fluid_model)
        s2 = Clapeyron.PH.entropy(fluid_model,p2,h2,z1)
        @test isapprox(s1,s2,atol=1e-5)
    end
end

@testset "Isentopic Expansion - Single Component" begin

    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]
        p1 = 101325*5;T1 = saturation_temperature(fluid_model,p1)[1] + 100.0; p2 = 101325;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= ThermoCycleGlides.isentropic_expander(p1,p2,1.0,h1,z1,fluid_model)
        s2 = Clapeyron.PH.entropy(fluid_model,p2,h2,z1)
        @test isapprox(s1,s2,atol=1e-5)
    end
end