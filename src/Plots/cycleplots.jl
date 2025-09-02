
"""
Plots TS diagram
"""
function plot_cycle(prob::HeatPump,sol::AbstractVector;N = 30,p_min = nothing)
  if isnothing(p_min)
      p_min = 0.5*sol[1]*101325
  end
    fig = plot()
    p_evap, p_cond = sol .* 101325 # convert to Pa
    fig_phase = plot_phase(fig,prob.fluid,prob.z; N = N, p_min = p_min)

    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_comp_in = enthalpy(prob.fluid, p_evap, T_evap_out,prob.z, phase = :vapor)
    h_comp_out = ThermoCycleGlides.isentropic_compressor(p_evap, p_cond, prob.η_comp, h_comp_in, prob.z, prob.fluid)
    h_comp_array = collect(range(h_comp_in, h_comp_out, length = N))
    p_comp_array = collect(range(p_evap, p_cond, length = N))
    T_ph(p,h) = Clapeyron.PH.temperature(prob.fluid, p, h, prob.z)
    T_comp_array = T_ph.(p_comp_array, h_comp_array)
    s_pt(p,t) = entropy(prob.fluid, p, t, prob.z, phase = :vapor)
    s_comp_array = s_pt.(p_comp_array, T_comp_array)
    plot!(fig_phase, s_comp_array, T_comp_array, label = "Compressor", color = :blue)

    T_cond_out = Clapeyron.bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_array = collect(range(h_cond_out, h_comp_out, length = N))
    T_cond_array = T_ph.(p_cond, h_cond_array)
    s_ph(p,h) = Clapeyron.PH.entropy(prob.fluid, p, h, prob.z)
    s_cond_array = s_ph.(p_cond, h_cond_array)
    plot!(fig_phase, s_cond_array, T_cond_array, label = "Condenser", color = :red)
    h_valve_in = h_cond_out
    h_valve_out = h_valve_in
    h_valve_array = collect(range(h_valve_in, h_valve_out, length = N))
    p_valve_array = collect(range(p_cond, p_evap, length = N))
    T_valve_array = T_ph.(p_valve_array, h_valve_array)
    s_valve_array = s_ph.(p_valve_array, h_valve_array)
    plot!(fig_phase, s_valve_array, T_valve_array, label = "Expansion Valve", color = :green)
    h_evap_in = h_valve_out
    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_evap_array = collect(range(h_evap_in, h_evap_out, length = N))
    T_evap_array = T_ph.(p_evap, h_evap_array)
    s_evap_array = s_ph.(p_evap, h_evap_array)
    plot!(fig_phase, s_evap_array, T_evap_array, label = "Evaporator", color = :purple)

    # secondary fluids
    T_cond_sf_array = collect(range(prob.T_cond_in, prob.T_cond_out, length = N))
    plot!(fig_phase,s_cond_array, T_cond_sf_array, label = "Condenser Secondary Fluid", color = :orange, linestyle = :dash)

    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length = N))
    plot!(fig_phase,s_evap_array, T_evap_sf_array, label = "Evaporator Secondary Fluid", color = :purple, linestyle = :dash)
    return fig_phase
end

function plot_cycle(prob::ORC, sol::AbstractVector;N = 30,p_min = nothing)
    if isnothing(p_min)
      p_min = 0.5*sol[2]*101325
    end
    fig = plot()
    p_evap, p_cond = sol .* 101325 # convert to Pa
    fig_phase = plot_phase(fig,prob; N = N, p_min = p_min)
    T_pump_in = bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_pump_in = Clapeyron.enthalpy(prob.fluid, p_cond, T_pump_in, prob.z)
    h_pump_out = ThermoCycleGlides.isentropic_pump(p_cond, p_evap, prob.η_pump, h_pump_in, prob.z, prob.fluid)
    h_pump_array = collect(range(h_pump_in, h_pump_out, length = N))
    T_ph(p,h) = Clapeyron.PH.temperature(prob.fluid, p, h, prob.z)
    T_pump_array = T_ph.(p_cond, h_pump_array)
    s_pt(p,t) = entropy(prob.fluid, p, t, prob.z)
    s_ph(p,h) = Clapeyron.PH.entropy(prob.fluid, p, h, prob.z)
    s_pump_array = s_pt.(p_cond, T_pump_array)
    plot!(fig_phase, s_pump_array, T_pump_array, label = "Pump", color = :blue)

    T_evap_out = dew_temperature(prob.fluid, p_evap, prob.z)[1] + prob.ΔT_sh
    h_evap_out = Clapeyron.enthalpy(prob.fluid, p_evap, T_evap_out, prob.z)
    h_evap_in = h_pump_out
    h_evap_array = collect(range(h_evap_in, h_evap_out, length = N))
    T_evap_array = T_ph.(p_evap, h_evap_array)
    s_evap_array = s_ph.(p_evap, h_evap_array)
    plot!(fig_phase, s_evap_array, T_evap_array, label = "Evaporator", color = :purple)

    h_exp_in = h_evap_out
    h_exp_out = ThermoCycleGlides.isentropic_expander(p_evap, p_cond, prob.η_expander, h_exp_in, prob.z, prob.fluid)
    h_exp_array = collect(range(h_exp_in, h_exp_out, length = N))
    p_exp_array = collect(range(p_evap, p_cond, length = N))
    T_exp_array = T_ph.(p_exp_array, h_exp_array)
    s_exp_array = s_ph.(p_exp_array, h_exp_array)
    plot!(fig_phase, s_exp_array, T_exp_array, label = "Expander", color = :orange)

    h_cond_in = h_exp_out
    T_cond_out = bubble_temperature(prob.fluid, p_cond, prob.z)[1] - prob.ΔT_sc
    h_cond_out = Clapeyron.enthalpy(prob.fluid, p_cond, T_cond_out, prob.z)
    h_cond_array = collect(range(h_cond_in, h_cond_out, length = N))
    T_cond_array = T_ph.(p_cond, h_cond_array)
    s_cond_array = s_ph.(p_cond, h_cond_array)
    plot!(fig_phase, s_cond_array, T_cond_array, label = "Condenser", color = :red)
    # secondary fluids
    T_evap_sf_array = collect(range(prob.T_evap_out, prob.T_evap_in, length = N))
    plot!(fig_phase, s_evap_array, T_evap_sf_array, label = "Evaporator Secondary Fluid", color = :purple, linestyle = :dash)
    T_cond_sf_array = collect(range(prob.T_cond_out, prob.T_cond_in, length = N))
    plot!(fig_phase, s_cond_array, T_cond_sf_array, label = "Condenser Secondary Fluid", color = :orange, linestyle = :dash)
    # @show minimum(T_evap_sf_array .- T_evap_array)
    return fig_phase
end

function plot_phase(fig::Plots.Plot,prob::ORC;N = 30,p_min = 101325*0.4)
  if length(prob.fluid.components) == 1
    return plot_phase_pure(fig,prob.fluid; N = N, p_min = p_min)
  end
  if length(prob.fluid.components) >1
    return plot_phase_mix(fig,prob.fluid,prob.z,N=N,p_min = p_min)
  end
  throw(error("IDK what you have passed to the function"))
end

function plot_phase(fig::Plots.Plot,prob::HeatPump; N = 30,p_min = 101325*0.4)

  if length(prob.fluid.components) == 1
    return plot_phase_pure(fig,prob.fluid; N = N, p_min = p_min)
  end
  if length(prob.fluid.components) >1
    return plot_phase_mix(fig,prob.fluid,prob.z,N=N,p_min = p_min)
  end
  throw(error("IDK what you have passed to the function"))
end

function plot_phase(fig::Plots.Plot,fluid::EoSModel,z::AbstractVector;N = 100,p_min = 101325*0.4)
  if length(z) == 1
    return plot_phase_pure(fig,fluid,N = N,p_min = p_min)
  end
  if length(z) > 1
    return plot_phase_mix(fig,fluid,z,N=N,p_min = p_min)
  end
  throw(error("IDK what you have passed to the function"))
end

function plot_phase_pure(fig::Plots.Plot,fluid::EoSModel; N = 100,p_min = 101325*0.4)
    @assert length(fluid.components) == 1
    pcrit = crit_pure(fluid)[2]
    p_array = collect(range(p_min, pcrit, length = N))
    Tsat(x) = saturation_temperature(fluid, x)[1]
    T = Tsat.(p_array)
    T[end] = crit_pure(fluid)[1]
    s_l = similar(T)
    for i in eachindex(T)
        s_l[i] = entropy(fluid,p_array[i],T[i],phase =:liquid)
    end
    s_g = similar(T)
    for i in eachindex(T)
        s_g[i] = entropy(fluid, p_array[i], T[i],phase =:vapor)
    end

    plot!(fig,s_l, T, label = false, ylabel = "Temperature (K)", xlabel = "Entropy (J/K)", title = "$(fluid.components[1])")
    plot!(fig, s_g, T, label = false)
    return fig
end


function plot_phase_mix(fig::Plots.Plot,fluid::EoSModel,z::AbstractVector;N = 100, p_min = 101325*0.4)
  @assert length(fluid.components) >= 2 
  Tcrit,pcrit,_ = crit_mix(fluid,z)
  p_array = collect(range(p_min, pcrit, length = N))
  Tdew(x) = dew_temperature(fluid,x,z)[1]
  Tbub(x) = bubble_temperature(fluid,x,z)[1]
  Td = Tdew.(p_array);
  Tb = Tbub.(p_array);
  Td[end] = Tcrit;
  Tb[end] = Tcrit;
  ThermoCycleGlides.fix_nan!(Td)
  ThermoCycleGlides.fix_nan!(Tb)
  sb = similar(Tb);
  for i in eachindex(Tb)
      sb[i] = entropy(fluid,p_array[i],Tb[i],z,phase =:liquid)
  end
  sd = similar(Td)
  for i in eachindex(Td)
    sd[i] = entropy(fluid, p_array[i], Td[i],z,phase =:vapor)
  end
  plot!(fig,sb, Tb, label = false, ylabel = "Temperature (K)", xlabel = "Entropy (J/K)", title = "$(fluid.components)")
  plot!(fig, sd, Td, label = false)
end

export plot_cycle, plot_phase