"""
optimize(prob::ThermoCycleProblem; kwargs...)

Its goal is to find the optimal subcooling and superheating values that maximize the cycle performance (COP or ORC -efficiency).

This should return the optimized cycle struct and the result from Metaheuristics
"""
function optimize(prob::ThermoCycleProblem,algo::Metaheuristics.Algorithm; autodiff::Bool = true,N::Int = 20,lower = 2.0,tol = 1e-3,xtol = 1e-6,ftol = 1e-6,max_iter= 1000)
        ub = max(abs(prob.T_cond_out - prob.T_evap_out),abs(prob.T_evap_out - prob.T_cond_in)) # upper bound on subcooling & superheating
        ub = ones(2) * ub
        lb = ones(2) * lower # lower bound on subcooling & superheating
        bounds  = Metaheuristics.BoxConstrainedSpace(lb = lb, ub = ub)
        soltemp = deepcopy(prob)
        function obj(x)
            try
                prob_ = deepcopy(prob)
                prob_.ΔT_sc = x[2]
                prob_.ΔT_sh = x[1]
                sol = solve(prob_,autodiff = autodiff,N = N,xtol = xtol,ftol = ftol,max_iter = max_iter)
                if norm(sol.residuals) > tol
                    return 0.0
                else
                    if prob isa HeatPump || prob isa HeatPumpRecuperator
                        if sol.x[1] > sol.x[2]
                            return 0.0
                        end
                        return COP(prob_,sol)
                    elseif prob isa ORC || prob isa ORCEconomizer
                        if sol.x[1] < sol.x[2]
                            return 0.0
                        end
                        return η(prob_,sol)
                    end
                end
            catch
                return 0.0
            end
        end
        if algo.options.parallel_evaluation == true
            function obj_parallel(x)
                fitness = zeros(size(X,1))
                Threads.@threads for i in eachindex(x)
                    fitness[i] = obj(x[i])
                end
                return fitness
            end
            result = Metaheuristics.optimize(obj_parallel,bounds,algo)
            soltemp.ΔT_sc = result.best_sol.x[2]
            soltemp.ΔT_sh = result.best_sol.x[1]
            return result,soltemp
        end

        if algo.options.parallel_evaluation == false
            result = Metaheuristics.optimize(obj,bounds,algo)
                        soltemp.ΔT_sc = result.best_sol.x[2]
            soltemp.ΔT_sh = result.best_sol.x[1]
            return result,soltemp
        end
end

function optimize(prob::ThermoCycleProblem,algo::Metaheuristics.Algorithm,params::ThermoCycleParameters)
    return optimize(prob,algo,N = params.N,autodiff = params.autodiff,xtol = params.xtol,ftol = params.ftol,max_iter = params.max_iters,tol = params.restart_TOL)
end

export optimize