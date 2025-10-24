# function F(prob::ThermoCycleProblem,x::AbstractVector{T},N::Int,internal_pinch::Bool) where T<:Real
#     return F(prob,x,N=N,internal_pinch = internal_pinch)
# end


function generate_initial_point(prob::HeatPump,lb::AbstractVector{T},ub::AbstractVector{T}) where T<: Real
    return [
        minimum(lb), maximum(ub)
    ]    
end

function generate_initial_point(prob::ORC,lb::AbstractVector{T},ub::AbstractVector{T}) where T<: Real
    return [
        maximum(ub), minimum(lb)
    ]    
end

function generate_initial_point(prob::ORCEconomizer,lb::AbstractVector{T},ub::AbstractVector{T}) where T<: Real
    return [
        maximum(ub), minimum(lb)
    ]    
end


function generate_initial_point(prob::HeatPumpRecuperator,lb::AbstractVector{T},ub::AbstractVector{T}) where T<: Real
    return [
        minimum(lb), maximum(ub)
    ]    
end

"""
`generate_box_solve_bounds(prob::HeatPump) -> lb, ub`
Generates lower and upper bounds for the heat pump problem based on its parameters.
"""
function generate_box_solve_bounds(prob::HeatPump)
    Tcrit,pcrit,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_cond_out > Tcrit
        psat_max = 0.99*pcrit 
    else 
        psat_max = bubble_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond + prob.ΔT_sc,prob.z)[1] 
        if !isfinite(psat_max)
            psat_max = 0.99*pcrit 
        end
    end

    psat_min = dew_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1]
    
    ub[1] = psat_max
    lb[1] = psat_min

    ub[2] = psat_max
    lb[2] = psat_min
    if !isfinite(psat_max)
        throw(error("The upper bound on pressure is not finite. Check cycle parameters. Possible error on bubble_pressure calculation."))
    end
    return lb./101325, ub./101325 # normalize to 101325 Pa
end

function generate_box_solve_bounds(prob::HeatPumpRecuperator)
    Tcrit,_,_ = crit_mix(prob.hp.fluid, prob.hp.z)
    lb = zeros(eltype(prob.hp.z), 2)
    ub = zeros(eltype(prob.hp.z), 2)
    if prob.hp.T_cond_out > Tcrit
        psat_max = 0.95*pcrit
    else
        psat_max = bubble_pressure(prob.hp.fluid,prob.hp.T_cond_out + prob.hp.pp_cond + prob.hp.ΔT_sc,prob.hp.z)[1]
        if !isfinite(psat_max)
            psat_max = 0.95*pcrit 
        end 
    end
    psat_min = dew_pressure(prob.hp.fluid,prob.hp.T_evap_out - prob.hp.pp_evap - prob.hp.ΔT_sh,prob.hp.z)[1]
    ub[1] = psat_max
    lb[1] = psat_min

    ub[2] = psat_max
    lb[2] = psat_min
    if !isfinite(psat_max)
        throw(error("The upper bound on pressure is not finite. Check cycle parameters. Possible error on bubble_pressure calculation."))
    end


    return lb./101325, ub./101325 # normalize to 101325 Pa
end

function generate_box_solve_bounds(prob::ORC)
    Tcrit,pcrit,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_evap_in > Tcrit
         psat_max = 0.95*pcrit
    end
    if prob.T_evap_out > Tcrit
        throw(error("For now we handel subcritical ORC. The outlet temperature of the evap : 
        $(prob.T_evap) is higher than critical temperature of the fluid $Tcrit, this will not allow to meet the pinch points.
        "))
    end
    psat_min = dew_pressure(prob.fluid,prob.T_cond_in,prob.z)[1]
    psat_max = bubble_pressure(prob.fluid,prob.T_evap_in - prob.ΔT_sh,prob.z)[1] #pcrit*0.9#dew_pressure(prob.fluid,prob.T_evap_in,prob.z)[1]
    ub[1] = psat_max#dew_pressure(prob.fluid,prob.T_evap_in - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    lb[1] = psat_min#bubble_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    ub[2] = psat_max
    lb[2] = psat_min 

    return lb./101325, ub./101325 # normalize to 101325 Pa
end

"""
Struct for solver and system paramters
    N -  hex discrtization (incase of mixtures)
    internal_pinch - to check if pinch point was met in two-phase
"""
mutable struct ThermoCycleParameters
    N::Int
    autodiff::Bool
    internal_pinch::Bool
    fd_order::Int
    xtol::Real
    ftol::Real
    restart_TOL::Real
    max_iters::Int
end

function ThermoCycleParameters(; 
    N::Int = 20,
    autodiff::Bool = true,
    internal_pinch::Bool = true,
    fd_order::Int = 2,
    xtol::Real = 1e-6,
    ftol::Real = 1e-6,
    restart_TOL::Real = 1e-3,
    max_iters::Int = 100
)
    N > 0 || error("N must be positive")
    if !autodiff && fd_order < 2
        error("If autodiff is false, fd_order must be ≥ 2 (higher-order finite differences required).")
    end
    return ThermoCycleParameters(N, autodiff, internal_pinch, fd_order, xtol, ftol, restart_TOL, max_iters)
end

# function check(prob::ThermoCycleProblem,param::ThermoCycleParameters)
#     if param.internal_pinch == true && length(prob.z) == 1
#         param.internal_pinch = false 
#     end
# end



function generate_box_solve_bounds(prob::ORCEconomizer)
    Tcrit,_,_ = crit_mix(prob.orc.fluid, prob.orc.z)
    lb = zeros(eltype(prob.orc.z), 2)
    ub = zeros(eltype(prob.orc.z), 2)
    if prob.orc.T_evap_in > Tcrit
        throw(error("For now only subcritical ORC are supported. Inlet temperature to evaporator must be below critical temperature."))
    end
    psat_min = dew_pressure(prob.orc.fluid,prob.orc.T_cond_in,prob.orc.z)[1]
    psat_max = bubble_pressure(prob.orc.fluid,prob.orc.T_evap_in,prob.orc.z)[1] #pcrit*0.9#dew_pressure(prob.fluid,prob.T_evap_in,prob.z)[1]
    ub[1] = psat_max#dew_pressure(prob.fluid,prob.T_evap_in - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    lb[1] = psat_min#bubble_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    ub[2] = psat_max
    lb[2] = psat_min 

    return lb./101325, ub./101325 # normalize to 101325 Pa
end

function solve_ad(prob::ThermoCycleProblem,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3,xtol = 1e-8,ftol = 1e-8,max_iter= 1000)
    f(x::AbstractVector{T}) where {T<:Real} = F(prob, x,N = N)
    T = promote_type(typeof(lb), typeof(ub))
    x0 = generate_initial_point(prob,lb,ub)
    sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter)
    sol.soltype = :subcritical
    if norm(sol.residuals)  > restart_TOL
        # f(x::AbstractVector{T}) where {T<:Real} = F_transcritical(prob, x,N = N)
        # x0 = generate_bounds(prob,lb,ub)
        # sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter)
        # sol.soltype = :transcritical
        x0 = (lb + ub) ./ 2
        sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter)
    end
    return sol
end

function solve_fd(prob::ThermoCycleProblem,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3,fd_order = 2,xtol = 1e-8,ftol = 1e-8,max_iter= 1000)
    f(x::AbstractVector{T}) where {T<:Real} = F(prob, x,N = N)
    x0 = generate_initial_point(prob,lb,ub)
    sol = constrained_newton_fd(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter,fd_order = fd_order)
    sol.soltype = :subcritical
    if norm(sol.residuals) > restart_TOL
        x0 = (lb + ub) ./ 2
        sol = constrained_newton_fd(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter,fd_order = fd_order)
    end
    return sol
end


"""
    Solves for pressure values in HP and ORC cycles for the given glide and problem parameters. 
    Define those problems in the respective structs. 
    For now the default box-nonlinear solver is newton-raphson, but this can be changed to other solvers in the future.
"""
function solve(prob::ThermoCycleProblem;autodiff::Bool = true, fd_order =2 , N::Int64 = 20,restart_TOL = 1e-3,xtol = 1e-6,ftol = 1e-6,max_iter= 1000)
    lb,ub = generate_box_solve_bounds(prob)
    if autodiff
        return sol = solve_ad(prob, lb, ub, N = N, restart_TOL = restart_TOL, xtol = xtol, ftol = ftol, max_iter = max_iter)
    else
        return sol = solve_fd(prob, lb, ub, N = N,fd_order = fd_order,restart_TOL = restart_TOL, xtol = xtol, ftol = ftol, max_iter = max_iter) 
    end
end

function solve(prob::ThermoCycleProblem,param::ThermoCycleParameters)
    return solve(prob,autodiff = param.autodiff,fd_order=param.fd_order,restart_TOL = param.restart_TOL,N = param.N,xtol = param.xtol,
    ftol = param.ftol,max_iter= param.max_iters)
end

export solve, ThermoCycleParameters

function show(io::IO,params::ThermoCycleParameters)
    println(io, "ThermoCycleParameters:")
    println(io, "  N               = ", params.N)
    println(io, "  autodiff        = ", params.autodiff)
    println(io, "  internal_pinch  = ", params.internal_pinch)
    println(io, "  fd_order        = ", params.fd_order)
    println(io, "  xtol            = ", params.xtol)
    println(io, "  ftol            = ", params.ftol)
    println(io, "  restart_TOL     = ", params.restart_TOL)
    println(io, "  max_iters       = ", params.max_iters)
end



