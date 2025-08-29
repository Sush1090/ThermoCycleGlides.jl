

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
    psat_max = bubble_pressure(prob.fluid,prob.T_cond_out + prob.pp_cond + prob.ΔT_sc,prob.z)[1] 
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

"""
returns the residues of the pinch points for `HeatPump`
"""
function F(prob::HeatPump,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    # @show x
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh 
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    h_comp_in = h_evap_out; 
    # @show h_comp_in, p_evap,p_cond
    h_comp_out = isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid)
    # @show h_comp_out
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_in = h_comp_out

    T_cond(h) = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)
    T_evap(h) = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)

    h_cond_array = collect(range(h_cond_out, h_cond_in, length=N))
    T_cond_array = T_cond.(h_cond_array)
    # fix_nan!(T_cond_array)
    T_cond_sf_array = collect(range(prob.T_cond_in, prob.T_cond_out, length=N))
    
    ΔTpp_cond = minimum(T_cond_array .- T_cond_sf_array) - prob.pp_cond
    h_valve_in = h_cond_out;
    h_valve_out = h_valve_in # isenthalpic expansion

    h_evap_array = collect(range(h_valve_out, h_evap_out, length=N))
    T_evap_array = T_evap.(h_evap_array)
    # fix_nan!(T_evap_array)
    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length=N))
    ΔTpp_evap = minimum(T_evap_sf_array .- T_evap_array) - prob.pp_evap

    return [ΔTpp_evap, ΔTpp_cond]
end

function F(prob::ORC,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_evap_out = Clapeyron.dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    # @show T_evap_out
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    h_exp_in = h_evap_out;
    h_exp_out = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)
    h_cond_in = h_exp_out
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_array = collect(range(h_cond_out, h_cond_in, length=N))
    T_cond(h) = Clapeyron.PH.temperature(prob.fluid, p_cond, h, prob.z)
    T_cond_array = T_cond.(h_cond_array)
    # fix_nan!(T_cond_array)
    T_cond_sf_array = collect(range(prob.T_cond_in, prob.T_cond_out, length=N))
    ΔTpp_cond = minimum(T_cond_array .- T_cond_sf_array) - prob.pp_cond
    h_pump_in = h_cond_out
    h_pump_out = ThermoCycleGlides.isentropic_compressor(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_evap_array = collect(range(h_pump_out, h_evap_out, length=N))
    T_evap(h) = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)
    T_evap_array = T_evap.(h_evap_array)
    # fix_nan!(T_evap_array)
    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length=N))
    ΔTpp_evap = minimum(T_evap_sf_array .- T_evap_array) - prob.pp_evap
    return [ΔTpp_evap, ΔTpp_cond]#, T_cond_array, T_evap_array, T_cond_sf_array, T_evap_sf_array
end

"""
Returns pressures and residues to the problem specificed HP/ORC using AD.
"""
function solve_ad(prob::HeatPump,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3)
    f(x::AbstractVector{T}) where {T<:Real} = F(prob, x,N = N)
    x0 =  zeros(2) #(lb + ub) ./ 2 # initial guess
    x0[1] = maximum(ub)
    x0[2] = minimum(lb)
    sol = constrained_newton_ad(f, x0, lb, ub; xtol = 1e-8, ftol = 1e-8, iterations = 100)
    if norm(f(sol)) > restart_TOL
        x0 = (lb + ub) ./ 2
        sol = constrained_newton_ad(f, x0, lb, ub; xtol = 1e-8, ftol = 1e-8, iterations = 100)
    end
    return sol, f(sol)
end


"""
Returns pressures and residues to the problem specificed HP/ORC using FD. 
"""
function solve_fd(prob::ThermoCycleGlides.HeatPump,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20)
    f(x::AbstractVector{T}) where {T<:Real} = ThermoCycleGlides.F(prob, x,N = N)
    x0 = (lb + ub) ./ 2 # initial guess
    sol = ThermoCycleGlides.constrained_newton_fd(f, x0, lb, ub; xtol = 1e-8, ftol = 1e-8, iterations = 100)
    return sol, f(sol)
end


function solve_ad(prob::ORC,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20,restart_TOL = 1e-3)
    f(x::AbstractVector{T}) where {T<:Real} = F(prob, x,N = N)
    x0 =  zeros(2) #(lb + ub) ./ 2 # initial guess
    x0[1] = maximum(ub)
    x0[2] = minimum(lb)
    sol = constrained_newton_ad(f, x0, lb, ub; xtol = 1e-8, ftol = 1e-8, iterations = 100)
    if norm(f(sol)) > restart_TOL
        x0 = (lb + ub) ./ 2
        sol = constrained_newton_ad(f, x0, lb, ub; xtol = 1e-8, ftol = 1e-8, iterations = 100)
    end
    return sol, f(sol)
end

function solve_fd(prob::ORC,lb::AbstractVector,ub::AbstractVector;N::Int64 = 20)
    f(x::AbstractVector{T}) where {T<:Real} = F(prob, x,N = N)
    # x0 = (lb + ub) ./ 2 # initial guess
    x0 =  zeros(2) #(lb + ub) ./ 2 # initial guess
    x0[1] = maximum(ub)
    x0[2] = minimum(lb)
    sol = constrained_newton_fd(f, x0, lb, ub; xtol = 1e-8, ftol = 1e-8, iterations = 100)
    return sol, f(sol)
end





"""
    Solves for pressure values in HP and ORC cycles for the given glide and problem parameters. 
    Define those problems in the respective structs. 
    For now the default box-nonlinear solver is newton-raphson, but this can be changed to other solvers in the future.
"""
function solve(prob::HeatPump;autodiff::Bool=true,N::Int64 = 20)
    lb,ub = generate_box_solve_bounds(prob)
    if autodiff
        return sol,res = solve_ad(prob, lb, ub, N = N)
    else
        return sol,res = solve_fd(prob, lb, ub, N = N)
    end
end

function solve(prob::ORC;autodiff::Bool=true,N::Int64 = 20)
    lb,ub = generate_box_solve_bounds(prob)
    if autodiff
        return sol,res = solve_ad(prob, lb, ub, N = N)
    else
        return sol,res = solve_fd(prob, lb, ub, N = N)
    end
end


export solve



