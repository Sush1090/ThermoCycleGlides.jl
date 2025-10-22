
"""
   ORC{T<:Real} <: ThermoCycleProblem

Defines an Organic Rankine Cycle (ORC) problem with thermodynamic and design 
parameters specified in Kelvin and dimensionless efficiencies.

# Fields
- `fluid::EoSModel`: Equation of State (EoS) model representing the working fluid.
- `z::AbstractVector{T}`: Mole fraction composition vector of the working fluid.
- `T_evap_in::T`: Inlet temperature of the evaporator [K].
- `T_evap_out::T`: Outlet temperature of the evaporator [K].
- `ΔT_sh::T`: Degree of superheating at the expander inlet [K].
- `T_cond_in::T`: Inlet temperature of the condenser [K].
- `T_cond_out::T`: Outlet temperature of the condenser [K].
- `ΔT_sc::T`: Degree of subcooling at the pump inlet [K].
- `η_pump::T`: Isentropic efficiency of the pump [-].
- `η_expander::T`: Isentropic efficiency of the expander [-].
- `pp_evap::T`: Minimum temperature difference (pinch point) at the evaporator [K].
- `pp_cond::T`: Minimum temperature difference (pinch point) at the condenser [K].

"""
mutable struct ORC{T<:Real} <: ThermoCycleProblem
    fluid::EoSModel
    z::AbstractVector{T}
    T_evap_in::T
    T_evap_out::T
    ΔT_sh::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc::T
    η_pump::T
    η_expander::T
    pp_evap::T
    pp_cond::T
end


function ORC(; fluid::EoSModel, z, T_evap_in, T_evap_out, T_cond_in, T_cond_out,
             η_pump, η_expander, pp_evap, pp_cond, ΔT_sh, ΔT_sc)

    @assert length(z) > 0 "Composition vector z must not be empty"
    @assert length(fluid.components) == length(z) "Composition vector z must match number of fluid components"

    # Heat source (secondary fluid) temperature drop
    @assert T_evap_in > T_evap_out "Evaporator secondary fluid must cool down (T_evap_in > T_evap_out)"

    # Heat sink (secondary fluid) temperature rise
    @assert T_cond_out > T_cond_in "Condenser secondary fluid must heat up (T_cond_out > T_cond_in)"

    # ORC thermodynamic requirement
    @assert T_evap_out > T_cond_in "Working fluid evaporation temperature must exceed condensation temperature"

    # Efficiency assertions
    @assert 0 < η_pump <= 1 "Pump efficiency must be in (0, 1]"
    @assert 0 < η_expander <= 1 "Expander efficiency must be in (0, 1]"

    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"

    @assert ΔT_sh ≥ 0 "Superheating must be non-negative"
    @assert ΔT_sc ≥ 0 "Subcooling must be non-negative"

    # Ensure subcritical evaporation
    # Tcrit, _, _ = crit_mix(fluid, z)
    # @assert T_evap_out < Tcrit - pp_evap "Evaporator outlet must be below critical temperature minus pinch"

    # Promote scalars
    type_promoted = promote_type(
        eltype(z),
        typeof(T_evap_in), typeof(T_evap_out),
        typeof(T_cond_in), typeof(T_cond_out),
        typeof(η_pump), typeof(η_expander),
        typeof(pp_evap), typeof(pp_cond),
        typeof(ΔT_sh), typeof(ΔT_sc)
    )

    z_T = convert(Vector{type_promoted}, z)
    T_evap_in_T  = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T  = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_pump_T     = convert(type_promoted, η_pump)
    η_expander_T = convert(type_promoted, η_expander)
    pp_evap_T    = convert(type_promoted, pp_evap)
    pp_cond_T    = convert(type_promoted, pp_cond)
    ΔT_sh_T      = convert(type_promoted, ΔT_sh)
    ΔT_sc_T      = convert(type_promoted, ΔT_sc)

    return ORC( 
        fluid,
        z_T,
        T_evap_in_T,
        T_evap_out_T,
        ΔT_sh_T,
        T_cond_in_T,
        T_cond_out_T,
        ΔT_sc_T,
        η_pump_T,
        η_expander_T,
        pp_evap_T,
        pp_cond_T
    )
end



"""
p[1] -> condensor pressure
p[2] -> evaporator pressure
"""
function η(prob::ORC,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    Tsat_cond = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1]
    h1 = Clapeyron.enthalpy(prob.fluid, p_cond, Tsat_cond - prob.ΔT_sc, prob.z)
    h2 = isentropic_pump(p_cond, p_evap, prob.η_pump, h1, prob.z, prob.fluid)

    Tsat_evap = Clapeyron.dew_temperature(prob.fluid, p_evap, prob.z)[1]
    h3 = Clapeyron.enthalpy(prob.fluid, p_evap, Tsat_evap + prob.ΔT_sh, prob.z)
    h4 = isentropic_expander(p_evap, p_cond, prob.η_expander, h3, prob.z, prob.fluid)
    return ((h4-h3) + (h2-h1))/(h3-h2)
end

"""
`show_parameters(prob::ORC)` prints parameters in REPL
"""
function show_parameters(prob::ORC)
    println("ORC Parameters:")
    println("Fluid: ", prob.fluid)
    println("Composition: ", prob.z)
    println("Evaporator Inlet Temperature: ", prob.T_evap_in, " K")
    println("Evaporator Outlet Temperature: ", prob.T_evap_out, " K")
    println("Condenser Inlet Temperature: ", prob.T_cond_in, " K")
    println("Condenser Outlet Temperature: ", prob.T_cond_out, " K")
    println("Pump Efficiency: ", prob.η_pump)
    println("Expander Efficiency: ", prob.η_expander)
    println("Evaporator Pinch Point: ", prob.pp_evap)
    println("Condenser Pinch Point: ", prob.pp_cond)
    println("Superheating: ", prob.ΔT_sh)
    println("Subcooling: ", prob.ΔT_sc)
end
export ORC, η


function F(prob::ORC,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    if length(prob.fluid.components) == 1
        return F_pure(prob,x)
    end
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
    h_pump_out = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_evap_array = collect(range(h_pump_out, h_evap_out, length=N))
    T_evap(h) = Clapeyron.PH.temperature(prob.fluid, p_evap, h, prob.z)
    T_evap_array = T_evap.(h_evap_array)
    # fix_nan!(T_evap_array)
    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length=N))
    ΔTpp_evap = minimum(T_evap_sf_array .- T_evap_array) - prob.pp_evap
    return [ΔTpp_evap, ΔTpp_cond]#, T_cond_array, T_evap_array, T_cond_sf_array, T_evap_sf_array
end


function F_pure(prob::ORC,x::AbstractVector{T}) where T<:Real
    @assert length(prob.fluid.components) == 1 "Pure fluid dispatch only"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_evap_out = Clapeyron.saturation_temperature(prob.fluid, p_evap)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_exp_in = h_evap_out;
    h_exp_out = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)
    h_cond_in = h_exp_out
    T_cond_out = Clapeyron.saturation_temperature(prob.fluid, p_cond)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    T_cond_sat = Clapeyron.saturation_temperature(prob.fluid, p_cond)[1]
    h_cond_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_sat,prob.z,phase = :liquid)
    h_cond_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_sat,prob.z,phase = :vapour)
    h_cond_array = [h_cond_in,h_cond_sat_vapour,h_cond_sat_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.fluid,p_cond,h_cond_array,prob.z)
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    
    h_pump_in = h_cond_out
    h_pump_out = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_evap_in = h_pump_out
    T_evap_sat = Clapeyron.saturation_temperature(prob.fluid, p_evap)[1]
    h_evap_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_sat,prob.z,phase = :liquid)
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_sat,prob.z,phase = :vapour)
    h_evap_array = [h_evap_out,h_evap_sat_vapour,h_evap_sat_liquid,h_pump_out]
    T_evap_array = Clapeyron.PH.temperature.(prob.fluid,p_evap,h_evap_array,prob.z)
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)

    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond

    return [ΔT_evap,ΔT_cond]
end

"""
Function that gives specific power ratings for ORC by fixing outlet power of expander to equal 1.
"""
function power_ratings(prob::ORC,sol::AbstractVector{T}) where T
    p_evap,p_cond = sol.*101325
    
    T_out_evap = dew_temperature(prob.fluid,p_evap,prob.z)[1] + prob.ΔT_sh
    h_out_evap = enthalpy(prob.fluid,p_evap,T_out_evap,prob.z)
    h_in_exp = h_out_evap;
    h_out_exp = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.η_expander, h_in_exp, prob.z, prob.fluid)
    Δh_exp = h_out_exp - h_in_exp
    if Δh_exp > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after expansion should be negative"
    end

    h_in_cond = h_out_exp
    T_out_cond = bubble_temperature(prob.fluid,p_cond,prob.z)[1] - prob.ΔT_sc
    h_out_cond = enthalpy(prob.fluid,p_cond,T_out_cond,prob.z)
    ΔQ_cond = h_out_cond - h_in_cond
    if ΔQ_cond > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after condensation should be negative"
    end
    h_in_pump = h_out_cond
    h_out_pump = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.η_pump, h_in_pump, prob.z, prob.fluid)
    Δh_pump = h_out_pump - h_in_pump
    if Δh_pump < 0 
         @warn "something wrong in the system. Change in enthalpy of the fluid after pump should be positive"
    end
    h_in_evap = h_out_pump
    ΔQ_evap = h_out_evap - h_in_evap
    if ΔQ_evap < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after evaporator should be positive"
    end

    return [Δh_exp/Δh_exp , Δh_pump/Δh_exp , ΔQ_evap/Δh_exp, ΔQ_cond/Δh_exp]
end

function power_ratings(prob::ORC,sol::SolutionState)
    return power_ratings(prob,sol.x)
end

"""
        ORCEconomizer{T<:Real} <: ThermoCycleProblem

Defines an Organic Rankine Cycle (ORC) configuration with an economiser (regenerative heat exchanger),
extending the base `ORC` problem with a specified effectiveness.

# Fields
- `orc::ORC{T}`: Base ORC system definition containing the thermodynamic parameters.
- `ϵ::T`: Effectiveness of the economiser (regenerator) [-].
"""
mutable struct ORCEconomizer{T<:Real} <: ThermoCycleProblem
    orc::ORC{T}
    ϵ::T
end

function ORCEconomizer(; orc::ORC, ϵ::Real)
    @assert 0.0 ≤ ϵ < 1.0 "Economizer effectiveness must be in [0, 1)"
    type_promoted = promote_type(eltype(orc.z), typeof(ϵ))
    ϵ_T = convert(type_promoted, ϵ)
    return ThermoCycleGlides.ORCEconomizer(orc, ϵ_T)
end

function F_pure(prob::ORCEconomizer,x::AbstractVector{T}) where T<:Real
    @assert length(prob.orc.fluid.components) == 1 "Pure fluid dispatch only"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_sat_evap = Clapeyron.saturation_temperature(prob.orc.fluid, p_evap)[1]
    T_sat_cond = Clapeyron.saturation_temperature(prob.orc.fluid, p_cond)[1]
    T_evap_out = Clapeyron.saturation_temperature(prob.orc.fluid, p_evap)[1] + prob.orc.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.orc.fluid, p_evap, T_evap_out, prob.orc.z)
    h_exp_in = h_evap_out;
    h_exp_out = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h_exp_in, prob.orc.z, prob.orc.fluid)
    T_exp_out = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h_exp_out, prob.orc.z)

    T_cond_out = Clapeyron.saturation_temperature(prob.orc.fluid, p_cond)[1] - prob.orc.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_cond_out, prob.orc.z)
    h_pump_in = h_cond_out
    h_pump_out = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h_pump_in, prob.orc.z, prob.orc.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h_pump_out, prob.orc.z)

    q_ihex = IHEX_Q(prob.orc.fluid, prob.ϵ,T_exp_out,p_cond,T_pump_out,p_evap,prob.orc.z)

    h_evap_in = h_pump_out + q_ihex
    h_cond_in = h_exp_out - q_ihex

    h_evap_sat_liquid = Clapeyron.enthalpy(prob.orc.fluid,p_evap,T_sat_evap,prob.orc.z,phase = :liquid)
    h_evap_sat_gas = Clapeyron.enthalpy(prob.orc.fluid,p_evap,T_sat_evap,prob.orc.z,phase = :vapour)
    h_evap_array = [h_evap_out,h_evap_sat_gas,h_evap_sat_liquid,h_evap_in]
    T_evap_array = Clapeyron.PH.temperature.(prob.orc.fluid,p_evap,h_evap_array,prob.orc.z)
    T_evap_sf_f(h) = prob.orc.T_evap_in - (h_evap_out - h)*(prob.orc.T_evap_in - prob.orc.T_evap_out)/(h_evap_out - h_evap_in)
    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.orc.pp_evap

    h_cond_sat_liquid = Clapeyron.enthalpy(prob.orc.fluid,p_cond,T_sat_cond,prob.orc.z,phase = :liquid)
    h_cond_sat_gas = Clapeyron.enthalpy(prob.orc.fluid,p_cond,T_sat_cond,prob.orc.z,phase = :vapour)
    h_cond_array = [h_cond_in,h_cond_sat_gas,h_cond_sat_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.orc.fluid,p_cond  ,h_cond_array,prob.orc.z)
    T_cond_sf_f(h) = prob.orc.T_cond_out - (h_cond_in - h)*(prob.orc.T_cond_out - prob.orc.T_cond_in)/(h_cond_in - h_cond_out)
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.orc.pp_cond
    return [ΔT_evap,ΔT_cond]
end



function F(prob::ORCEconomizer,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    if length(prob.orc.fluid.components) == 1
        return F_pure(prob,x)
    end
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_sat_evap = Clapeyron.dew_temperature(prob.orc.fluid, p_evap, prob.orc.z)[1]
    T_sat_cond = Clapeyron.bubble_temperature(prob.orc.fluid, p_cond, prob.orc.z)[1]
    T_evap_out = T_sat_evap + prob.orc.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.orc.fluid, p_evap, T_evap_out, prob.orc.z)
    h_exp_in = h_evap_out;
    h_exp_out = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h_exp_in, prob.orc.z, prob.orc.fluid)
    T_exp_out = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h_exp_out, prob.orc.z)
    T_cond_out = T_sat_cond - prob.orc.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.orc.fluid, p_cond, T_cond_out, prob.orc.z)
    h_pump_in = h_cond_out
    h_pump_out = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h_pump_in, prob.orc.z, prob.orc.fluid)
    T_pump_out = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h_pump_out, prob.orc.z)
    q_ihex = IHEX_Q(prob.orc.fluid, prob.ϵ,T_exp_out,p_cond,T_pump_out,p_evap,prob.orc.z)
    h_evap_in = h_pump_out + q_ihex
    h_cond_in = h_exp_out - q_ihex
    h_evap_array = collect(range(h_evap_out, h_evap_in, length=N))
    T_evap(h) = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h, prob.orc.z)
    T_evap_array = T_evap.(h_evap_array)
    T_evap_sf(h) = prob.orc.T_evap_in - (h_evap_out - h)*(prob.orc.T_evap_in - prob.orc.T_evap_out)/(h_evap_out - h_evap_in)
    T_evap_sf_array = T_evap_sf.(h_evap_array)
    ΔTpp_evap = minimum(T_evap_sf_array .- T_evap_array) - prob.orc.pp_evap
    h_cond_array = collect(range(h_cond_out, h_cond_in, length=N))
    T_cond(h) = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h, prob.orc.z)
    T_cond_array = T_cond.(h_cond_array)
    T_cond_sf(h) = prob.orc.T_cond_in + (h - h_cond_out)*(prob.orc.T_cond_out - prob.orc.T_cond_in)/(h_cond_in - h_cond_out)
    T_cond_sf_array = T_cond_sf.(h_cond_array)
    ΔTpp_cond = minimum(T_cond_array .- T_cond_sf_array) - prob.orc.pp_cond
    return [ΔTpp_evap, ΔTpp_cond]#, T_cond_array,
end


function η(prob::ThermoCycleGlides.ORCEconomizer,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    Tsat_cond = Clapeyron.bubble_temperature(prob.orc.fluid, p_cond, prob.orc.z)[1]
    h1 = Clapeyron.enthalpy(prob.orc.fluid, p_cond, Tsat_cond - prob.orc.ΔT_sc, prob.orc.z)
    h2 = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.orc.η_pump, h1, prob.orc.z, prob.orc.fluid)
    T2 = Clapeyron.PH.temperature(prob.orc.fluid, p_evap, h2, prob.orc.z)
    Tsat_evap = Clapeyron.dew_temperature(prob.orc.fluid, p_evap, prob.orc.z)[1]
    h3 = Clapeyron.enthalpy(prob.orc.fluid, p_evap, Tsat_evap + prob.orc.ΔT_sh, prob.orc.z)
    h4 = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.orc.η_expander, h3, prob.orc.z, prob.orc.fluid)
    T4 = Clapeyron.PH.temperature(prob.orc.fluid, p_cond, h4, prob.orc.z)
    q_ihex = ThermoCycleGlides.IHEX_Q(prob.orc.fluid, prob.ϵ,T4,p_cond,T2,p_evap,prob.orc.z)
    h2_new = h2 + q_ihex
    return ((h4-h3) + (h2-h1))/(h3-h2_new)
end

export ORCEconomizer


"""
    η(prob::ThermoCycleGlides.ThermoCycleProblem, sol::SolutionState) -> Float64

Computes the thermal efficiency of a thermodynamic cycle given a problem definition 
and its corresponding solution state.

# Arguments
- `prob::ThermoCycleGlides.ThermoCycleProblem`: The thermodynamic cycle problem 
  containing fluid properties, boundary conditions, and component parameters.
- `sol::SolutionState`: The solution state object containing the converged 
  state variables (`x`), residuals, and convergence information.

# Returns
- `Float64`: The computed cycle efficiency, defined as the ratio of net work 
  output to heat input.

# Notes
This method acts as a wrapper that extracts the solution vector `x` from 
`sol` and calls the lower-level `η(prob, x)` implementation.
"""
function η(prob::ThermoCycleGlides.ThermoCycleProblem,sol::SolutionState)
    return η(prob,sol.x)
end



function _F(prob::ORC, x::AbstractVector{T}) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    @assert length(prob.z) != 1 "This implementation is for mixtures"

    p_evap,p_cond = x .* 101325 # convert to Pa
    T_evap_out = Clapeyron.dew_temperature(prob.fluid, p_evap,prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_exp_in = h_evap_out;
    h_exp_out = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)
    h_cond_in = h_exp_out
    T_cond_in = Clapeyron.PH.temperature(prob.fluid,p_cond,h_cond_in,prob.z)
    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond,prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    T_cond_dew = Clapeyron.dew_temperature(prob.fluid, p_cond,prob.z)[1]
    T_cond_bub = Clapeyron.bubble_temperature(prob.fluid, p_cond,prob.z)[1]
    h_cond_sat_liquid = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_bub,prob.z,phase = :liquid)
    h_cond_sat_vapour = Clapeyron.enthalpy(prob.fluid,p_cond,T_cond_dew,prob.z,phase = :vapour)
    h_cond_array = [h_cond_in,h_cond_sat_vapour,h_cond_sat_liquid,h_cond_out]
    T_cond_array = [T_cond_in,T_cond_dew,T_cond_bub,T_cond_out]
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    
    h_pump_in = h_cond_out
    h_pump_out = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_evap_in = h_pump_out
    T_evap_in = Clapeyron.PH.temperature(prob.fluid,p_evap,h_evap_in,prob.z)[1]
    T_evap_bub = Clapeyron.bubble_temperature(prob.fluid, p_evap,prob.z)[1]
    T_evap_dew = Clapeyron.dew_temperature(prob.fluid, p_evap,prob.z)[1]
    T_evap_out = T_evap_dew + prob.ΔT_sh
    h_evap_bub = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_bub,prob.z)
    h_evap_dew = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_dew,prob.z)
    h_evap_out = Clapeyron.enthalpy(prob.fluid,p_evap,T_evap_out,prob.z)
    h_evap_array = [h_evap_out,h_evap_bub,h_evap_dew,h_evap_in]
    T_evap_array = zeros(length(h_evap_array))
    for i in eachindex(h_evap_array)
        T_evap_array[i] = Clapeyron.PH.temperature(prob.fluid,p_evap,h_evap_array[i],prob.z)
    end
    
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)

    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond
    return [ΔT_evap,ΔT_cond]
end