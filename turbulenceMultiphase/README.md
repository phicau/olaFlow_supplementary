olaFlow_supplementary - turbulenceMultiphase
======

# Description

OpenFOAM® does not provide by default incompressible turbulence models for multiphase systems (i.e. models that take into account the density variation between the air and water phases). This often results in an excessive wave damping as the simulation progresses, due to turbulence build-up (increasing *nut*).

This repository provides modified versions of k-ε and k-ω SST turbulence models to simulate correctly multiphase systems and mitigate this effect (*kEpsilonMultiphase*, *kOmegaSSTMultiphase*).

Moreover, a new variation of the k-ω SST model was recently proposed in an excellent paper by Devolder et al. (2017) to model free surface cases (e.g. waves) correctly. The new model (*kOmegaSSTBuoyancy*) features an additional buoyancy term that helps suppress the spurious turbulence generation at the free surface.

The present implementation is (at least) compatible with OpenFOAM 5.0.0 and OpenFOAM v1712.

All the turbulence models in this library have been developed from the default turbulence models available in OpenFOAM and are solely based on the references given below.

Although the authors of the paper have released their turbulence model in GitHub, which can be found here: https://github.com/BrechtDevolder-UGent-KULeuven/buoyancyModifiedTurbulenceModels, this implementation has not been used in the development of the present library.

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

# Contents

## turbulenceMultiphaseLibrary

The library includes an *allMake* script for automatic compilation.

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

- How to create custom turbulence models: http://hassankassem.me/posts/newturbulencemodel/

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.