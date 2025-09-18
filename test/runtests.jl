using ThermoCycleGlides, Clapeyron, LinearAlgebra
using Test


include("fluids.jl")

# @testset "ThermoCycleGlides.jl" begin
#     # Write your tests here.
# end


@testset "Isentopic Compression - Single Component" begin

    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]; p1 = 101325;
        T1 = dew_temperature(fluid_model,p1,z1)[1] + 20;  p2 = 101325*2;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= ThermoCycleGlides.isentropic_compressor(p1,p2,1.0,h1,z1,fluid_model)
        s2 = Clapeyron.PH.entropy(fluid_model,p2,h2,z1)
        @test isapprox(s1,s2,atol=1e-5)
    end
end

@testset "Isentopic pump - Single Component" begin

    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]; p1 = 101325;
        T1 = bubble_temperature(fluid_model,p1,z1)[1] - 20;  p2 = 101325*2;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= ThermoCycleGlides.isentropic_pump(p1,p2,1.0,h1,z1,fluid_model)
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

@testset "ORC - mixture Non-linear solver" begin
    fluid = cPR(["propane","butane"],idealmodel = ReidIdeal);
    T_evap_out = 310
    z = [0.6157894736842106, 9.384210526315789]
    ΔT_sh = 11.473684210526315 
    _orc_ = ORC(fluid=fluid, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
    sol,res = solve(_orc_,N = 30)
    @test norm(res) < 1e-3
end



@testset "Pure - fluids - easy ORC NL solver" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 310
        ΔT_sh = 5.0 
        _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol,res = solve(_orc_,N = 30)
        @test norm(res) < 1e-3
    end
end


@testset "ORC - mixture Non-linear solver FD" begin
    fluid = cPR(["propane","butane"],idealmodel = ReidIdeal);
    T_evap_out = 310
    z = [0.6157894736842106, 9.384210526315789]
    ΔT_sh = 11.473684210526315 
    _orc_ = ORC(fluid=fluid, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
    sol,res = solve(_orc_,N = 30,autodiff = false)
    @test norm(res) < 1e-3
end



@testset "Pure - fluids - easy ORC NL solver" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 310
        ΔT_sh = 5.0 
        _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol,res = solve(_orc_,N = 30,autodiff = true)
        @test norm(res) < 1e-3
    end
end

@testset "Pure - fluids - hard ORC NL solver" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 340
        ΔT_sh = 5.0 
        _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=360, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol,res = solve(_orc_,N = 30,autodiff = false)
        @test norm(res) < 1e-3
    end
end

@testset "Pure - fluids - hard ORC-Economizer NL solver" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 340
        ΔT_sh = 5.0 
        _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=360, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol,res = solve(_orc_,N = 30,autodiff = false)
        ϵ = 0.7
        _orc_econ_ = ORCEconomizer(_orc_,ϵ)
        sol_e,res_e = solve(_orc_econ_,N = 30,autodiff = false)
        
        @test norm(res_e) < 1e-3
        @test abs(η(_orc_econ_,sol_e)) >= abs(η(_orc_,sol))
    end
end

@testset "Pure - fluids - easy HP NL solver" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 300
        ΔT_sh = 5.0 
        _hp_ = HeatPump(fluid=fluid_model, z=z, T_evap_in=310, T_evap_out=T_evap_out, T_cond_in=340, T_cond_out=355, η_comp=0.75, pp_evap=2, pp_cond=2, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol,res = solve(_hp_,N = 30,autodiff = false)
        @test norm(res) < 1e-3
    end
end


@testset "Pure : HeatPumpRecuperator solver" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 300
        ΔT_sh = 5.0 
        hp_ = HeatPump(fluid=fluid_model, z=z, T_evap_in=310, T_evap_out=T_evap_out, T_cond_in=340, T_cond_out=355, η_comp=0.75, pp_evap=2, pp_cond=2, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        ϵ = 0.4
        sol_hp,res_hp = solve(hp_,N = 30,autodiff = true)
        _hp_ = HeatPumpRecuperator(hp_,ϵ)
        sol,res = solve(_hp_,N = 30,autodiff = true)
        @test norm(res) < 1e-3
        # @test abs(COP(_hp_,sol)) >= abs(COP(hp_,sol_hp))
    end
end