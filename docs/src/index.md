
# ThermoCycleGlides.jl

Documentation of ThermoCycleGlides.jl.

The goal of this package is to have the non-linear solver for Heat Pump and Organic Rankine Cycle systems. It solves for pressures at evaporator and condensor for given pinch-point temperatures as well as provides plotting framework of the solution. 
As of now the package is robust for subcritical cycle parameters. For thermodynamic properties, [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl) is used as backend. 

For details of modeling see:
[Carnot batteries for heat and power coupling: Energy, Exergy, Economic and Environmental (4E) analysis - Laterre, Antoine](https://hdl.handle.net/2268/333259), which is describes modeling for pure fluids while this package also implements it for mixtures. 

The nonlinear solver chosen is Newton-Raphson with box bounds which is inspired by the implementation in NLboxsolve.jl

This is for steady-state applications.

There are 4 systems provided: 

1. Organic Rankine Cycle : `ORC`
2. Organic Rankine Cycle with internal heat exchanger : `ORCEconomizer`
3. Heat Pump : `HeatPump`
4. Heat Pump with internal heat exchanger : `HeatPumpRecuperator` 

The implemented version of these systems consist of the following components:
    
1. Compressor : Modeled with isentropic efficiency   
2. Expander : Modeled with isentropic efficiency
3. Value : Modeled as isenthalpic process
4. Evaporator : Volumes of equal change in enthalpy.
5. Condensor : Volumes of equal change in enthalpy.
6. Heat Exchangers (no phase change) : using $\epsilon$ as effectivness of heat exchanger. 




