/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "IncompressibleMomentumTransportModel.H"
#if OFVERSION >= 900
    #include "kinematicTransportModel.H"
#else 
    #include "transportModel.H"
#endif
#include "addToRunTimeSelectionTable.H"
#include "makeMomentumTransportModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeMomentumTransportModelTypes(Alpha, Rho, baseModel, BaseModel, Transport) \
namespace Foam                                                                 \
{                                                                              \
        typedef BaseModel<Transport> Transport##BaseModel;                     \
        typedef laminarModel<Transport##BaseModel>                             \
            laminar##Transport##BaseModel;                                     \
        typedef RASModel<Transport##BaseModel> RAS##Transport##BaseModel;      \
        typedef LESModel<Transport##BaseModel> LES##Transport##BaseModel;      \
}

makeMomentumTransportModelTypes
(
    geometricOneField,
    geometricOneField,
    incompressibleMomentumTransportModel,
    IncompressibleMomentumTransportModel,
    #if OFVERSION >= 900
        kinematicTransportModel
    #else 
        transportModel
    #endif
);

#if OFVERSION >= 900
    #define makeRASModel(Type)                                                 \
    makeTemplatedMomentumTransportModel                                        \
    (kinematicTransportModelIncompressibleMomentumTransportModel, RAS, Type)
#else
    #define makeRASModel(Type)                                                 \
    makeTemplatedMomentumTransportModel                                        \
    (transportModelIncompressibleMomentumTransportModel, RAS, Type)
#endif

// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "kEpsilonMultiphase.H"
makeRASModel(kEpsilonMultiphase);

#include "kOmegaSSTMultiphase.H"
makeRASModel(kOmegaSSTMultiphase);

#include "kOmegaSSTBuoyancy.H"
makeRASModel(kOmegaSSTBuoyancy);

#include "kOmegaSSTStable.H"
makeRASModel(kOmegaSSTStable);

// ************************************************************************* //
