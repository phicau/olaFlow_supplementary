olaFlow_supplementary - turbulenceMultiphase
======

# Description

OpenFOAM® does not provide by default incompressible turbulence models for multiphase systems (i.e. models that take into account the density variation between the air and water phases). Moreover, two-equation turbulence models have recently been proven unstable by Larsen & Fuhrman (2018). This often results in an excessive wave damping as the simulation progresses, due to turbulence build-up (increasing *nut*).

This repository provides modified versions of k-ε and k-ω SST turbulence models to simulate correctly multiphase systems, accounting for the density of each phase in the equations (*kEpsilonMultiphase*, *kOmegaSSTMultiphase*).

Moreover, new variations of the k-ω SST model were recently proposed in two excellent papers by Devolder et al. (2017, 2018) and Larsen & Fuhrman (2018) to model free surface cases (e.g. waves) correctly. The new models (*kOmegaSSTBuoyancy*, *kOmegaSSTStable*) feature an additional buoyancy term that helps suppress the spurious turbulence generation at the free surface. *kOmegaSSTStable* also includes a modified *nut* formulation to convert mitigate the instability.

The present implementation is (at least) compatible with OpenFOAM 8/dev and OpenFOAM 2006+.

All the turbulence models in this library have been developed from the default turbulence models available in OpenFOAM, based on the references given.

Devolder et al. have released their turbulence model in GitHub in which they also include a modified k-ω model and can be found here: https://github.com/BrechtDevolder-UGent-KULeuven/buoyancyModifiedTurbulenceModels  

Larsen & Fuhrman have released their turbulence model in GitHub in which they also include other modified models and can be found here: 
https://github.com/BjarkeEltardLarsen/RANS_stableOF50

# Usage

To use the new turbulence models in your cases, you need to include the following code in *controlDict*:
```
libs
(
    "libturbulenceMultiphaseOlaFlowModels.so"
);
```

Then, select either of the following in *turbulenceProperties*:

- kEpsilonMultiphase
- kOmegaSSTMultiphase
- kOmegaSSTBuoyancy
- kOmegaSSTStable

# Contents

## turbulenceMultiphaseLibrary_*

The library includes an *allMake* script for automatic compilation.  
Currently compatible with all OpenFOAM versions up to 8 and 2006+.

## Tutorials/turbulenceMultiphaseFlume

Regular waves propagating on a flume using *olaFlow*; very small turbulence levels should be expected.
To simulate the case run any of the *runCase...* scripts. Each one will simulate the case with a different turbulence model, as referenced in their names.

Videos of the comparison between the different turbulence models can be accessed at:
- https://youtu.be/20cpgpBtcGE
- https://youtu.be/eHLF_wzAMn4
- https://youtu.be/TPbodvlpLI8

# References

- **Application of a buoyancy-modified k-ω SST turbulence model to simulate wave run-up around a monopile subjected to regular waves using OpenFOAM**  
Brecht Devolder, Pieter Rauwoens & Peter Troch  
*Coastal Engineering* (2017), vol. 125, pp. 81–94 (https://doi.org/10.1016/j.coastaleng.2017.04.004)

- **Performance of a buoyancy-modified k-ω and k-ω SST turbulence model for simulating wave breaking under regular waves using OpenFOAM®**  
Brecht Devolder, Peter Troch & Pieter Rauwoens  
*Coastal Engineering* (2018), vol. 138, pp. 49–65 (https://doi.org/10.1016/j.coastaleng.2018.04.011)

- **On the over-production of turbulence beneath surface waves in RANS models**  
Larsen, B.E. & Fuhrman, D.R.  
*Journal of Fluid Mechanics* (2018), vol. 853, pp. 419-460
(https://doi.org/10.1017/jfm.2018.577)

- How to create custom turbulence models: http://hassankassem.me/posts/newturbulencemodel/

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.
