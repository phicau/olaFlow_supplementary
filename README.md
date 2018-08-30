olaFlow_supplementary
======

# Description

Supplementary materials for the olaFlow project (https://github.com/phicau/olaFlow), including new solvers, modules and tutorials.

All tutorials include *run* and *clean* scripts for automatic running/cleaning the case.

Compatibility has been tested with the OpenFOAM® versions reported only.

# Download and compilation - Basic download guide

To get a full copy of the olaFlow supplementary materials, run in a terminal:

`git clone git://github.com/phicau/olaFlow_supplementary.git`

Updates can be downloaded in the future from the *olaFlow_supplementary* folder as follows:

`git checkout`  
`git pull`

# Cases

## multiphaseInterFoam

- Multi-fluid tutorials (air-water-oil) with olaFlow waves.
- Cases to be run with OpenFOAM 5.0 vanilla solver *multiphaseInterFoam*.

## olaIsoFlow

- Coupling of *isoAdvector* (OpenFOAM 5.0 and v1706-v1712) sharp interface advection method with olaFlow wave boundary conditions.
- The tutorial case is a basic wave flume in 2D.
- More information can be found in: https://sites.google.com/view/olaflowcfd/blog/olaflow-coupling-with-isoadvector

## turbulenceMultiphase

- Library that implements density-aware turbulence models for multiphase modelling and a stable version of the k-ω SST model. These developments prevent turbulence buid-up to a large extent in wave simulations with olaFlow/OpenFOAM.
- The tutorial case is a basic wave flume in 2D.
- For further details regarding the implementation see the readme file in the turbulenceMultiphase folder.
- More information can be found in: https://sites.google.com/view/olaflowcfd/blog/turbulence-models-for-wave-simulations

The new models are:

- kEpsilonMultiphase
- kOmegaSSTMultiphase
- kOmegaSSTBuoyancy
- kOmegaSSTStable

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.