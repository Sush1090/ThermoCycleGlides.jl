
using PrecompileTools: @setup_workload, @compile_workload  


    @setup_workload begin
        @compile_workload begin
            fluid = cPR("propane",idealmodel = ReidIdeal)
            hp = HeatPump(fluid=fluid, z=[1.0], T_evap_in=310, T_evap_out=300.0, T_cond_in=325, T_cond_out=355, η_comp=0.75, pp_evap=5, pp_cond=5, ΔT_sc=5.0, ΔT_sh=5.0)
            sol = solve(hp,autodiff = true,N = 20,xtol = 1e-10,ftol = 1e-10,max_iter= 1000);
            sol_fd = solve(hp,autodiff = false,N = 20,xtol = 1e-10,ftol = 1e-10,max_iter= 1000);
        end

        @compile_workload  begin
            fluid = cPR(["propane","butane"],idealmodel = ReidIdeal)
            hp = HeatPump(fluid=fluid, z=[1.0,1.0], T_evap_in=310, T_evap_out=300.0, T_cond_in=325, T_cond_out=355, η_comp=0.75, pp_evap=5, pp_cond=5, ΔT_sc=5.0, ΔT_sh=5.0)
            sol = solve(hp,autodiff = true,N = 20,xtol = 1e-10,ftol = 1e-10,max_iter= 1000);
            sol_fd = solve(hp,autodiff = false,N = 20,xtol = 1e-10,ftol = 1e-10,max_iter= 1000);
        end

        @compile_workload  begin
            fluid_model = cPR("propane",idealmodel = ReidIdeal)
            z = [1.0]
            T_evap_out = 340
            ΔT_sh = 5.0 
            _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=360, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
            sol = solve(_orc_,N = 30,autodiff = false)
            ϵ = 0.7
            _orc_econ_ = ORCEconomizer(_orc_,ϵ)
            sol_e = solve(_orc_econ_,N = 30,autodiff = false,ftol = 1e-8,xtol = 1e-8)
            sol = solve(_orc_,N = 30,autodiff = true,ftol = 1e-8,xtol = 1e-8)
        end
                @compile_workload  begin
            fluid_model = cPR(["propane","butane"],idealmodel = ReidIdeal)
            z = [1.0,1.0]
            T_evap_out = 340
            ΔT_sh = 5.0 
            _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=360, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
            sol = solve(_orc_,N = 30,autodiff = false)
            ϵ = 0.7
            _orc_econ_ = ORCEconomizer(_orc_,ϵ)
            sol_e = solve(_orc_econ_,N = 30,autodiff = false,ftol = 1e-8,xtol = 1e-8)
            sol = solve(_orc_,N = 30,autodiff = true,ftol = 1e-8,xtol = 1e-8)
        end
    end

