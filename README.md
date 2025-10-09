# ThermoCycleGlides.jl

[![Build Status](https://github.com/Sush1090/ThermoCycleGlides.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/ThermoCycleGlides.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package aims to solve Heat Pump and ORC systems for given known temperature glides. For now it is robust for sub-critical cycles. 

The thermodynamic computations use Clapeyron.jl. 

# Usage
Usage Heat Pump Example :

```julia
julia> using Clapeyron, ThermoCycleGlides

julia> fluid = cPR(["cyclopentane"],idealmodel = ReidIdeal);

julia> η_comp = 0.75; pp_cond = 2; pp_evap = 2;

julia> T_evap_in = 273.15 + 10; T_evap_out = 273.15 + 0; T_cond_in = 273.15 + 50;  T_cond_out = 273.15+60;

julia> ΔT_sc = 3; ΔT_sh = 10;

julia> hp = HeatPump(fluid=fluid,z=[1.0],T_evap_in=T_evap_in,T_evap_out = T_evap_out,T_cond_in = T_cond_in,T_cond_out=T_cond_out,η_comp=η_comp,pp_evap=pp_evap,pp_cond=pp_cond,ΔT_sc = ΔT_sc,ΔT_sh = ΔT_sh);

julia> sol_hp = solve(hp)
SolutionState{Float64, Int64}([0.12829257763187535, 1.4551588056895837], 16, 7, [0.0, 5.684341886080802e-14], [0.07660159441435545, 0.07660159441435545], [1.6566058479359296, 1.6566058479359296], true, 0, 2.5468671054250572e-15, 8.038873388460929e-14)
 
julia> COP(hp,sol_hp)
-3.735868783526992
```

To plot do the following;

```julia
plot_cycle(hp,sol_hp,N=300)
```

![HP_cyclopentane](Images/hp_cyclopentane.png)


ORC Example:

```julia
julia> fluid = cPR(["propane"],idealmodel = ReidIdeal);

julia> orc = ORC(fluid = fluid,z = [1.0], T_evap_in = 360, T_evap_out = 340, T_cond_in = 280, T_cond_out = 290, η_expander = 0.75, η_pump = 0.8, ΔT_sh = 7.0, ΔT_sc= 3.0, pp_evap = 3.0, pp_cond = 3)
ORC{Float64}(PR{ReidIdeal, TwuAlpha, NoTranslation, vdW1fRule}("Propane"), [1.0], 360.0, 340.0, 7.0, 280.0, 290.0, 3.0, 0.8, 0.75, 3.0, 3.0)

julia> sol = solve(orc)
SolutionState{Float64, Int64}([27.24650730298343, 8.310317977639695], 10, 4, [5.684341886080802e-14, -1.1368683772161603e-13], [5.8848010540583715, 5.8848010540583715], [31.19314765995761, 31.19314765995761], true, 0, 8.321091160516346e-7, 3.7257249981036476e-6)

julia> η(orc,sol)
-0.08885630484821327
```

To plot the ORC cycle: 

```julia
julia> plot_cycle(orc,sol,N=300)
```

![orc_propane](Images/orc_propane.png)

To see the details of `SolutionState` do the following

```julia
julia> show_parameters(sol)
Iterations: 4
Function calls: 10
Final residual norm: 1.2710574864626038e-13
Final x: [27.24650730298343, 8.310317977639695]
Final lenx: 8.321091160516346e-7
Final lenf: 3.7257249981036476e-6
Lower bounds: [5.8848010540583715, 5.8848010540583715]
Upper bounds: [31.19314765995761, 31.19314765995761]
Autodiff: true
```