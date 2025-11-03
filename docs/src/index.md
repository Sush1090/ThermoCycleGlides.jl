
# ThermoCycleGlides.jl

Documentation of ThermoCycleGlides.jl.

The goal of this package is to have the non-linear solver for Heat Pump and Organic Rankine Cycle systems. This is for steady-state applications.
The simplest version of these systems consist of:
    
1. Compressor : Modeled with isentropic efficiency   
2. Expander : Modeled with isentropic efficiency
3. Value : Modeled as isenthalpic process
4. Evaporator : Volumes of equal change in enthalpy.
5. Condensor : Volumes of equal change in enthalpy.
6. Heat Exchangers (no phase change) : using $\epsilon$ as effectivness of heat exchanger. 




