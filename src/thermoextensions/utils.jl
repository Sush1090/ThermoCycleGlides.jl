

function isentropic_compressor(p_in::T1, p_out::T2, η_isen::T3, h_in::T4, z::AbstractArray, fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in, z)
    h_isen = Clapeyron.PS.enthalpy(fluid, p_out, s_isen, z)
    ha =  h_in + (h_isen - h_in) / η_isen
     T_out = Clapeyron.PH.temperature(fluid,p_out,ha,z)
     T_dew = dew_temperature(fluid,p_out,z)[1]
     if T_out < T_dew
        @warn "Fixing outlet of compressor at saturation temperature"
        return enthalpy(fluid,p_out,T_dew,z)
    end
    return ha
end

function isentropic_pump(p_in::T1, p_out::T2, η_isen::T3, h_in::T4, z::AbstractArray, fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in, z)
    h_isen = Clapeyron.PS.enthalpy(fluid, p_out, s_isen, z)
    return h_in + (h_isen - h_in) / η_isen
end

function isentropic_expander(p_in::T1,p_out::T2,η_isen::T3,h_in::T4,z::AbstractVector,fluid::EoSModel) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    s_isen = Clapeyron.PH.entropy(fluid, p_in, h_in,z) 
    h_isen = Clapeyron.PS.enthalpy(fluid, p_out, s_isen,z)
    h_out = h_in - (h_in - h_isen) * η_isen
    # force outlet to be gaseous 
    dt = Clapeyron.dew_temperature(fluid, p_out, z)[1]
    h_gas = Clapeyron.enthalpy(fluid, p_out, dt, z,phase = :gas)
    # @show   values(dt), h_gas, h_out, p_in, p_out
    if h_out < h_gas
        # @warn "The outlet enthalpy is below the gas enthalpy. Adjusting to gas phase."
        h_out = h_gas 
    end

    return h_out
end

function IHEX_Q(fluid::EoSModel,ϵ::T,T_in_left::T,p_in_left::T,T_in_right::T,p_in_right::T,z::AbstractVector{T}) where T<:Real
    @assert T_in_left >= T_in_right "Function assumse left fluid is hotter"
    c_in_left = Clapeyron.isobaric_heat_capacity(fluid,p_in_left,T_in_left,z)
    c_in_right = Clapeyron.isobaric_heat_capacity(fluid,p_in_right,T_in_right,z)
    C = min(c_in_left,c_in_right)
    Qmax = C*(T_in_left - T_in_right)
    Q = Qmax*ϵ
    return Q
end