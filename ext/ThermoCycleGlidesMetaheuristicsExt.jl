module CarnotMetaheuristicsExt

using Carnot
using Metaheuristics
using Carnot.CommonSolve


function _build_residual(prob::HeatPump, N::Int)
    f = let prob = prob, N = N
        x -> begin 
            nc = length(prob.fluid.components)
            return nc == 1 ? F_pure(prob, x) : F(prob, x, N = N)
        end
    end
end

function objective(prob::HeatPump,param::ThermoCycleParameters,x::AbstractVector)
    @assert length(x) == 2 "Only super and sub cool temperatures"
    prob.ΔT_sh = x[1]
    prob.ΔT_sc = x[2]
    show(param)
    #try 
    sol =  solve(prob,param)
        if Carnot.norm(sol.residuals) > 1e-3
            return 0.0
        else
            @show COP(prob,sol)
            return COP(prob,sol)
        end
    #catch 
      #  return 0.0
   # end
end

function _build_objective(prob::HeatPump,param::ThermoCycleParameters)
    f = let prob = prob
        x -> begin
            objective(prob,param,x)
        end 
    end
end
export _build_objective
function generate_optimization_bounds(prob::HeatPump)
    ΔT_sh_min = 0.0
    ΔT_sh_max = prob.T_evap_in - prob.T_evap_out
    ΔT_sc_min = 0.0
    ΔT_sc_max = prob.T_cond_out - prob.T_cond_in
    lb = [ΔT_sh_min,ΔT_sc_min]
    ub = [ΔT_sh_max,ΔT_sc_max]
    return lb,ub
end

function CommonSolve.solve(prob::HeatPump,
    alg::Metaheuristics.AbstractAlgorithm,param::ThermoCycleParameters)

    @time "Building Objective function..." begin
    ℓ = _build_objective(prob,param)
    end
    @time "Generating bounds ..." begin
    lb,ub = generate_optimization_bounds(prob)
    end
    # x01 = Carnot.generate_initial_point(prob,lb,ub,:average)
    # x02 = Carnot.generate_initial_point(prob,lb,ub,:default)
    # Metaheuristics.set_user_solutions!(alg, x01, ℓ)
    # Metaheuristics.set_user_solutions!(alg, x02, ℓ)
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    @info "Starting Optimization..."
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)
    
    x_best = Metaheuristics.minimizer(opt_result)

    loss_opt_M = Metaheuristics.minimum(opt_result)
  #  status_code = Symbol(Metaheuristics.TerminationStatusCode(opt_result))
    f_calls = Metaheuristics.nfes(opt_result)
    residuals = NaN .* x_best #don't calculate residuals
    result = SolutionState(
        x_best,                 #x::Vector{T}
        f_calls,                #f_calls::I
        opt_result.iteration,   #iterations::I
        residuals,              #residuals::Vector{T}
        lb,                     #lb::Vector{T}
        ub,                     #ub::Vector{T}
        false,                  #autodiff::Bool
        0,                      #fd_order::I
        0.0,                    #lenx::T
        loss_opt_M,             #lenf::T
        :optimized_solution             #soltype::Symbol
    )
end


# function CommonSolve.solve(prob::Carnot.HeatPump,
#     alg::Metaheuristics.AbstractAlgorithm,param::ThermoCycleParameters)

    
#     ℓ = Carnot._build_least_squares_objective(prob,N)
#     lb,ub = Carnot.generate_box_solve_bounds(prob)
#     x01 = Carnot.generate_initial_point(prob,lb,ub,:average)
#     x02 = Carnot.generate_initial_point(prob,lb,ub,:default)
#     Metaheuristics.set_user_solutions!(alg, x01, ℓ)
#     Metaheuristics.set_user_solutions!(alg, x02, ℓ)
#     bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
#     opt_result = Metaheuristics.optimize(ℓ,bounds,alg)
    
#     x_best = Metaheuristics.minimizer(opt_result)

#     loss_opt_M = Metaheuristics.minimum(opt_result)
#     status_code = Symbol(Metaheuristics.TerminationStatusCode(opt_result))
#     f_calls = Metaheuristics.nfes(opt_result)
#     residuals = NaN .* x_best #don't calculate residuals
#     result = SolutionState(
#         x_best,                 #x::Vector{T}
#         f_calls,                #f_calls::I
#         opt_result.iteration,   #iterations::I
#         residuals,              #residuals::Vector{T}
#         lb,                     #lb::Vector{T}
#         ub,                     #ub::Vector{T}
#         false,                  #autodiff::Bool
#         0,                      #fd_order::I
#         0.0,                    #lenx::T
#         loss_opt_M,             #lenf::T
#         status_code             #soltype::Symbol
#     )
# end

end #module