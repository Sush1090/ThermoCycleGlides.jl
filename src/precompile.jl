
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
    end

