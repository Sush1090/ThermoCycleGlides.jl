
mutable struct ORC{T<:Real}
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
    h2 = isentropic_compressor(p_cond, p_evap, prob.η_pump, h1, prob.z, prob.fluid)

    Tsat_evap = Clapeyron.dew_temperature(prob.fluid, p_evap, prob.z)[1]
    h3 = Clapeyron.enthalpy(prob.fluid, p_evap, Tsat_evap + prob.ΔT_sh, prob.z)
    h4 = isentropic_expander(p_evap, p_cond, prob.η_expander, h3, prob.z, prob.fluid)
    return ((h4-h3) - (h2-h1))/(h3-h2)
end


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
