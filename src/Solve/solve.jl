
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
    Tcrit,_,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_cond_out > Tcrit
        throw(error("For now only subcritical heat pumps are supported. Outlet temperature to condenser must be below critical temperature."))
    end
    psat_min = dew_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1]
    psat_max = bubble_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond,prob.z)[1] 
    ub[1] = psat_max#dew_pressure(prob.fluid,prob.T_evap_in - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    lb[1] = psat_min#bubble_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure

    ub[2] = psat_max#dew_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond,prob.z)[1] # condensor pressure
    lb[2] = psat_min #bubble_pressure(prob.fluid,prob.T_cond_in + prob.pp_cond + prob.ΔT_sc,prob.z)[1] # condensor pressure
    return lb./101325, ub./101325 # normalize to 101325 Pa
end

function generate_box_solve_bounds(prob::HeatPumpRecuperator)
    Tcrit,_,_ = crit_mix(prob.hp.fluid, prob.hp.z)
    lb = zeros(eltype(prob.hp.z), 2)
    ub = zeros(eltype(prob.hp.z), 2)
    if prob.hp.T_cond_out > Tcrit
        throw(error("For now only subcritical heat pumps are supported. Outlet temperature to condenser must be below critical temperature."))
    end
    psat_min = dew_pressure(prob.hp.fluid,prob.hp.T_evap_out - prob.hp.pp_evap - prob.hp.ΔT_sh,prob.hp.z)[1]
    psat_max = bubble_pressure(prob.hp.fluid,prob.hp.T_cond_out + prob.hp.pp_cond,prob.hp.z)[1]
    ub[1] = psat_max#dew_pressure(prob.fluid,prob.T_evap_in - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure
    lb[1] = psat_min#bubble_pressure(prob.fluid,prob.T_evap_out - prob.pp_evap - prob.ΔT_sh,prob.z)[1] # evaporator pressure

    ub[2] = psat_max#dew_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond,prob.z)[1] # condensor pressure
    lb[2] = psat_min #bubble_pressure(prob.fluid,prob.T_cond_in + prob.pp_cond + prob.ΔT_sc,prob.z)[1] # condensor pressure
    return lb./101325, ub./101325 # normalize to 101325 Pa
end

function generate_box_solve_bounds(prob::ORC)
    Tcrit,_,_ = crit_mix(prob.fluid, prob.z)
    lb = zeros(eltype(prob.z), 2)
    ub = zeros(eltype(prob.z), 2)
    if prob.T_evap_in > Tcrit
        throw(error("For now only subcritical ORC are supported. Inlet temperature to evaporator must be below critical temperature."))
    end
    psat_min = dew_pressure(prob.fluid,prob.T_cond_in,prob.z)[1]
    psat_max = bubble_pressure(prob.fluid,prob.T_evap_in,prob.z)[1] #pcrit*0.9#dew_pressure(prob.fluid,prob.T_evap_in,prob.z)[1]
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
    if norm(f(sol)) > restart_TOL
        x0 = (lb + ub) ./ 2
        sol = constrained_newton_ad(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter)
    end
    return sol, f(sol)
end

function solve_fd(prob::ThermoCycleProblem,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3,fd_order = 2,xtol = 1e-8,ftol = 1e-8,max_iter= 1000)
    f(x::AbstractVector{T}) where {T<:Real} = F(prob, x,N = N)
    x0 = generate_initial_point(prob,lb,ub)
    sol = constrained_newton_fd(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter,fd_order = fd_order)
    if norm(f(sol)) > restart_TOL
        x0 = (lb + ub) ./ 2
        sol = constrained_newton_fd(f, x0, lb, ub; xtol = xtol, ftol = ftol, iterations = max_iter,fd_order = fd_order)
    end
    return sol, f(sol)
end



"""
    Solves for pressure values in HP and ORC cycles for the given glide and problem parameters. 
    Define those problems in the respective structs. 
    For now the default box-nonlinear solver is newton-raphson, but this can be changed to other solvers in the future.
"""
function solve(prob::ThermoCycleProblem;autodiff::Bool = true, fd_order =2 , N::Int64 = 20,restart_TOL = 1e-3,xtol = 1e-8,ftol = 1e-8,max_iter= 1000)
        lb,ub = generate_box_solve_bounds(prob)
    if autodiff
        return sol,res = solve_ad(prob, lb, ub, N = N, restart_TOL = restart_TOL, xtol = xtol, ftol = ftol, max_iter = max_iter)
    else
        return sol,res = solve_fd(prob, lb, ub, N = N,fd_order = fd_order,restart_TOL = restart_TOL, xtol = xtol, ftol = ftol, max_iter = max_iter) 
    end
end

export solve



