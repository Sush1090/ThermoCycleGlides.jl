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
This package models and solves simple Heat pump (HP) and Organic Rankine Cycle (ORC) systems.  It solves the system given the secondary fluid temperature profiles and component parameters. The thermodynamic properties are computed using Clapeyron.jl[@Clapeyron-2022]. This allows modelling of this system with zeotropic mixtures, where the phase-change is non-isothermal unlike for pure fluids. The package implements a nonlinear pinch-point solver that determines the evaporator and condenser pressures corresponding to prescribed inlet and outlet temperatures, ensuring thermodynamic feasibility and consistent glide matching.  


# Statement of need
Designing and optimizing high-level HP-ORC systems requires handeling of temperature glides in the heat exchangers. This is generally done by  pinch-point analysis which finds the operating pressures in evaporators and condensors. This is a common technique that allows for mapping the performance of the cycles (for example used in: [] [] []). 

    

# Example
Should I put? Docs should be sufficient? Or just an image?
 

# References