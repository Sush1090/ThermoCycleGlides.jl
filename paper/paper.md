---
title: 'ThermoCycleGlides.jl: A package for solving simple HP-ORC systems'
tags:
  - Julia
  - Themodynamic
  - Heat Pump
  - ORC
  - Pinch Point
authors:
  - name: Sushrut Deshpande
    orcid: 0009-0004-6960-7641
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - given-names: Bart Janssens
    orcid: 0000-0001-6164-7893
    affiliation: 1
affiliations:
 - name: Royal Military Academy, Belgium
   index: 1
   ror: 02vmnye06
 - name: University of Liege, Belgium
   index: 2
date: 05 November 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary
This package models and solves simple Heat pump (HP) and Organic Rankine Cycle (ORC) systems.  It solves the system given the secondary fluid temperature profiles and component parameters. The thermodynamic properties are computed using [@Clapeyron-2022]. This allows modelling of this system with zeotropic mixtures, where the phase-change is non-isothermal unlike for pure fluids. The package implements a nonlinear pinch-point solver that determines the evaporator and condenser pressures corresponding to prescribed inlet and outlet temperatures, ensuring thermodynamic feasibility and consistent glide matching.  


# Statement of need
Thermodynamics is a vast discipline with applications that span scales—from quantum phenomena to large-scale energy systems. In practical energy engineering, thermodynamic principles underpin systems that rely on the behaviour of fluids to convert or transfer energy. Such systems are found in power generation and thermal management, where heat engines convert thermal energy into mechanical work, and heat pumps or refrigeration systems perform the inverse process, converting electrical energy into useful heating or cooling.

Both applications operate as thermodynamic cycles, in which a working fluid circulates through a series of components—such as compressors, expanders, evaporators, and condensers—to accomplish a specific objective. In a heat pump, the goal is to achieve a temperature lift, transferring heat from a lower-temperature source to a higher-temperature sink by supplying external energy. Conversely, in a heat engine, the aim is to extract work by transferring energy from a high-temperature source to a low-temperature sink.


Designing and optimising high-level heat pump (HP) and organic Rankine cycle (ORC) systems requires accurate handling of the secondary fluid temperature glides within the heat exchangers. These glides define how the external heat-source and heat-sink temperatures vary along the heat exchanger length and must be matched consistently with the working fluid’s thermodynamic behaviour. This is typically achieved through pinch-point analysis, which determines the evaporator and condenser operating pressures that satisfy thermal feasibility and ensure proper temperature matching.

ThermoCycleGlides.jl provides a computational framework for solving this nonlinear pinch-point problem and determining the corresponding system conditions. It also includes a plotting framework for visualising the thermodynamic cycle and glide matching. The package follows the system characterisation outlined by [@Antoine-Laterre-2025] and employs a discretised heat exchanger model assuming equal heat transfer in each control volume, as described in [@ZUHLSDORF].

    
# Example
I will change this part
![hp_example](./images/HP_example.png) 

![hp_example_sol](./images/hp_cyclopentane.png) 

# References