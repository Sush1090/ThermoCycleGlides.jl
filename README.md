# ThermoCycleGlides

[![Build Status](https://github.com/Sush1090/ThermoCycleGlides.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/ThermoCycleGlides.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package aims to solve Heat Pump and ORC systems for given known temperature glides. 

The thermodynamic computations use Clapeyron.jl. 

Usage Heat Pump Example :

```julia
julia> using Clapeyron, ThermoCycleGlides

julia> fluid = cPR(["cyclopentane"],idealmodel = ReidIdeal);

julia> η_comp = 0.75; pp_cond = 2; pp_evap = 2;

julia> T_evap_in = 273.15 + 10; T_evap_out = 273.15 + 0; T_cond_in = 273.15 + 50;  T_cond_out = 273.15+60;

julia> ΔT_sc = 3; ΔT_sh = 10;

julia> hp = HeatPump(fluid=fluid,z=[1.0],T_evap_in=T_evap_in,T_evap_out = T_evap_out,T_cond_in = T_cond_in,T_cond_out=T_cond_out,η_comp=η_comp,pp_evap=pp_evap,pp_cond=pp_cond,ΔT_sc = ΔT_sc,ΔT_sh = ΔT_sh);

julia> sol_hp, res_hp = solve(hp) # res_hp is the pinchpoint temperture residue.
([0.12829257763187563, 1.4440205301792532], [0.0, -4.6397769835948566e-5])
 
julia> COP(hp,sol_hp) # Computes COP of Heat pump
-3.7495668569405063
```

To plot do the following;

```julia
plot_cycle(hp,sol_hp,p_min=0.3*sol_hp[1]*101325,N = 1000)
```

![HP_cyclopentane](Images/hp_cyclopentane.png)


ORC Example:

```julia
julia> orc = ORC(fluid = fluid,z = [1.0], T_evap_in = 360, T_evap_out = 340, T_cond_in = 280, T_cond_out = 290, η_expander = 0.75, η_pump = 0.8, ΔT_sh = 7.0, ΔT_sc= 3.0, pp_evap = 3.0, pp_cond = 3)
ORC{Float64}(PR{ReidIdeal, TwuAlpha, NoTranslation, vdW1fRule}("Propane"), [1.0], 360.0, 340.0, 7.0, 280.0, 290.0, 3.0, 0.8, 0.75, 3.0, 3.0)

julia> sol,res = solve(orc)
([27.619094241454576, 8.297312696256173], [1.3659347928296484e-6, -3.1898690622256254e-6])

julia> η(orc,sol)
-0.11240385922458294
```

To plot the ORC cycle: 

```julia
julia> plot_cycle(orc,sol,p_min=0.5*sol_hp[2]*101325,N = 1000)
```

![orc_propane](Images/orc_propane.png)