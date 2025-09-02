
mutable struct HeatPump{T<:Real} <: ThermoCycleProblem
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


mutable struct HeatPumpRecuperator{T<:Real}
    hp::HeatPump{T}
    ϵ::T
end

function HeatPumpRecuperator(;hp,ϵ)
    @assert ϵ <= 1 "Effictiveness has to be less than 1"
    @assert ϵ >=0 "Effictiveness has to be greater than 0"
    return HeatPumpRecuperator(hp,ϵ)
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



"""
returns the residues of the pinch points for `HeatPump`
"""
function F(prob::HeatPump,x::AbstractVector{T};N::Int64) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    if length(prob.fluid.components) == 1
        return F_pure(prob,x)
    end
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

function F_pure(prob::HeatPump,x::AbstractVector{T}) where {T<:Real}
    @assert length(x) == 2 "x must be a vector of length 2"
    p_evap,p_cond = x .* 101325 # convert to Pa
    T_sat_evap = saturation_temperature(prob.fluid,p_evap)[1]
    T_sat_cond = saturation_temperature(prob.fluid,p_cond)[1]
    T_evap_out = T_sat_evap + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_comp_in = h_evap_out;
    h_comp_out = ThermoCycleGlides.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid)
    T_cond_out = T_sat_cond - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_in = h_comp_out
    h_cond_vapour = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_cond, prob.z,phase =:vapour)
    h_cond_liquid = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_cond, prob.z,phase =:liquid)
    h_cond_array = [h_cond_in,h_cond_vapour,h_cond_liquid,h_cond_out]
    T_cond_array = Clapeyron.PH.temperature.(prob.fluid,p_cond,h_cond_array,prob.z)
    T_cond_sf_f(h) = prob.T_cond_out - (h_cond_in - h)*(prob.T_cond_out - prob.T_cond_in)/(h_cond_in - h_cond_out)
    ΔT_cond = minimum(T_cond_array .- T_cond_sf_f.(h_cond_array)) - prob.pp_cond
    h_valve_in = h_cond_out;
    h_valve_out = h_valve_in # isenthalpic expansion

    h_evap_in = h_valve_out
    h_evap_sat_vapour = Clapeyron.enthalpy(prob.fluid, p_cond, T_sat_evap, prob.z,phase =:vapour)
    h_evap_array = reverse([h_evap_in,h_evap_sat_vapour,h_evap_out])
    T_evap_array = Clapeyron.PH.temperature.(prob.fluid,p_evap,h_evap_array,prob.z)
    T_evap_sf_f(h) = prob.T_evap_in - (h_evap_out - h)*(prob.T_evap_in - prob.T_evap_out)/(h_evap_out - h_evap_in)
    ΔT_evap = minimum(T_evap_sf_f.(h_evap_array) .- T_evap_array) - prob.pp_evap
    [ΔT_cond,ΔT_evap]
end

"""
Function that gives specific power ratings for HP by fixing outlet power of compressor to equal 1.
"""
function power_ratings(prob::HeatPump,sol::AbstractVector{T}) where T
    p_evap,p_cond = sol.*101325
    
    T_out_evap = dew_temperature(prob.fluid,p_evap,prob.z)[1] + prob.ΔT_sh
    h_out_evap = enthalpy(prob.fluid,p_evap,T_out_evap,prob.z)
    h_in_comp = h_out_evap;
    h_out_comp = ThermoCycleGlides.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_in_comp, prob.z, prob.fluid)
    Δh_comp = h_out_comp - h_in_comp
    if Δh_comp < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after compressor should be positive"
    end

    h_in_cond = h_out_comp
    T_out_cond = bubble_temperature(prob.fluid,p_cond,prob.z)[1] - prob.ΔT_sc
    h_out_cond = enthalpy(prob.fluid,p_cond,T_out_cond,prob.z)
    ΔQ_cond = h_out_cond - h_in_cond
    if ΔQ_cond > 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after condensation should be negative"
    end
    h_in_valve = h_out_cond
    h_out_valve = h_in_valve
    Δh_valve = h_in_valve - h_out_valve

    h_in_evap = h_out_valve
    ΔQ_evap = h_out_evap - h_in_evap
    if ΔQ_evap < 0
        @warn "something wrong in the system. Change in enthalpy of the fluid after evaporator should be positive"
    end

    return [Δh_comp/Δh_comp , Δh_valve/Δh_comp , ΔQ_evap/Δh_comp, ΔQ_cond/Δh_comp]
end

export power_ratings