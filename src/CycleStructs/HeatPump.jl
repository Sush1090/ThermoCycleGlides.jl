
mutable struct HeatPump{T<:Real}
    fluid::EoSModel
    z::AbstractVector{T}
    T_evap_in::T
    T_evap_out::T
    ΔT_sh::T
    T_cond_in::T
    T_cond_out::T
    ΔT_sc::T
    η_comp::T
    pp_evap::T
    pp_cond::T
end

function HeatPump(;fluid::EoSModel,z,T_evap_in,T_evap_out,T_cond_in,T_cond_out,η_comp,pp_evap,pp_cond,ΔT_sh,ΔT_sc)

    #default assertions
    @assert length(z) > 0 "Composition vector z must not be empty"
    @assert length(fluid.components) == length(z) "Composition vector z must match the number of components in the fluid model"
    @assert T_evap_in > T_evap_out "Evaporator inlet temperature must be more than outlet temperature for the secondary fluid"
    @assert T_cond_in < T_cond_out "Condenser inlet temperature must be less than outlet temperature for the secondary fluid"
    @assert η_comp > 0 && η_comp <= 1 "Compressor efficiency must be between 0 and 1"
    @assert pp_evap > 0 "Evaporator pinch point must be positive"
    @assert pp_cond > 0 "Condenser pinch point must be positive"
    @assert ΔT_sh >= 0 "Superheating temperature must be non-negative"
    @assert ΔT_sc >= 0 "Subcooling temperature must be non-negative"

    # Thermodynamic assertions
    # For heat-pump the inlet temperature of the condensor should be higher than inelt temperature of the evaporator
    @assert T_cond_in > T_evap_out "Condenser inlet temperature must be higher than evaporator outlet temperature for the heat pump to function properly"
    # inlet temperature of the condensor should be subcritical - pinch point
    Tcrit,_,_ = crit_mix(fluid,z)
    @assert T_cond_in < Tcrit - pp_cond "Condenser inlet temperature must be less than critical temperature ($Tcrit) minus pinch point ($pp_cond) for the heat pump to function properly"



    type_promoted = promote_type(eltype(z), typeof(T_evap_in), typeof(T_evap_out), typeof(T_cond_in), typeof(T_cond_out), typeof(η_comp), typeof(pp_evap), typeof(pp_cond), typeof(ΔT_sh), typeof(ΔT_sc))
    z_T = convert(Vector{type_promoted}, z)
    ΔT_sh_T = convert(type_promoted, ΔT_sh)
    ΔT_sc_T = convert(type_promoted, ΔT_sc)
    T_evap_in_T = convert(type_promoted, T_evap_in)
    T_evap_out_T = convert(type_promoted, T_evap_out)
    T_cond_in_T = convert(type_promoted, T_cond_in)
    T_cond_out_T = convert(type_promoted, T_cond_out)
    η_comp_T = convert(type_promoted, η_comp)
    pp_evap_T = convert(type_promoted, pp_evap)
    pp_cond_T = convert(type_promoted, pp_cond)
    return HeatPump(
    fluid,         # EoSModel
    z_T,             # AbstractVector{T}
    T_evap_in_T,   # T
    T_evap_out_T,  # T
    ΔT_sh_T,       # T
    T_cond_in_T,   # T
    T_cond_out_T,  # T
    ΔT_sc_T,       # T
    η_comp_T,      # T
    pp_evap_T,     # T
    pp_cond_T      # T
)
end


function COP(prob::HeatPump,sol::AbstractVector{T}) where {T<:Real}
    @assert length(sol) == 2 "Pressure vector p must be of length 2"
    p_evap,p_cond = sol .* 101325 # convert to Pa
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh 
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)

    h_comp_in = h_evap_out; 
    h_comp_out = isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid)

    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)

    return  (h_cond_out - h_comp_out)/(h_comp_out - h_comp_in) 
end


mutable struct HeatPumpRecuperator{T}
    hp::HeatPump{T}
    ϵ::T
end

function show_parameters(prob::HeatPump)
    println("Heat Pump Parameters:")
    println("Fluid: ", prob.fluid)
    println("Composition: ", prob.z)
    println("Evaporator Inlet Temperature: ", prob.T_evap_in, " K")
    println("Evaporator Outlet Temperature: ", prob.T_evap_out, " K")
    println("Superheating Temperature: ", prob.ΔT_sh, " K")
    println("Condenser Inlet Temperature: ", prob.T_cond_in, " K")
    println("Condenser Outlet Temperature: ", prob.T_cond_out, " K")
    println("Subcooling Temperature: ", prob.ΔT_sc, " K")
    println("Compressor Efficiency: ", prob.η_comp)
    println("Evaporator Pinch Point: ", prob.pp_evap)
    println("Condenser Pinch Point: ", prob.pp_cond)
end

export HeatPump, HeatPumpRecuperator, COP, show_parameters
