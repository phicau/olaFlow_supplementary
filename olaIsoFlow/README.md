olaFlow_supplementary - olaIsoFlow
======

# Description

- Coupling of *isoAdvector* (OpenFOAM 5.0 and v1706-v1712) sharp interface advection method with olaFlow wave boundary conditions.
- The tutorial case is a basic wave flume in 2D.
- More information can be found in: https://sites.google.com/view/olaflowcfd/blog/olaflow-coupling-with-isoadvector

# Contents

## solver_OF5

*olaIsoFlow* solver source code for OpenFOAM 5.0.
You are required to download and compile *isoAdvector*, as instructed in https://github.com/isoAdvector/isoAdvector
To compile the solver, run the *allMake* script. You will be prompted to input the location of the *isoAdvector* source code folder (e.g. /home/user/isoAdvector/OpenFOAM-5.0/src).

## solver_OFv17xx

*olaIsoFlow* solver source code for OpenFOAM v1706-v1712.
To compile run the *allMake* script.

## Tutorials/waveFlume_isoAdvector

Regular waves on a flume using isoAdvector. To simulate the case run the *runCase* script.
The video of the completed case can be accessed at: https://youtu.be/AjSYfB63LpI

----------------------------------------------------------
OPENFOAM®  is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com.