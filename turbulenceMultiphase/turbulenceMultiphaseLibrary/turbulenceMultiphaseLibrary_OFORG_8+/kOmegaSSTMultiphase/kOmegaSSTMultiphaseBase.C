/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "kOmegaSSTMultiphaseBase.H"
#if OFVERSION >= 900
    #include "fvModels.H"
    #include "fvConstraints.H"
#else 
    #include "fvOptions.H"
#endif
#include "bound.H"
#include "wallDist.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSSTMultiphase::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(muTP()/rhoTP())/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSSTMultiphase::
F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(muTP()/rhoTP())/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSSTMultiphase::
F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(muTP()/rhoTP())/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSSTMultiphase::
F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    #if OFVERSION >= 900
        fvConstraints::New(this->mesh_).constrain(this->nut_);
    #else
        fv::options::New(this->mesh_).correct(this->nut_);
    #endif
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::
correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*rhoTP()().dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::
omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*rhoTP()().dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*rhoTP()().dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSSTMultiphase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    MomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
bool kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::read()
{
    if (MomentumTransportModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSSTMultiphase<MomentumTransportModel, BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    
    // Local references
    const alphaField& alpha = this->alpha_;
    // const rhoField& rho = this->rho_;
    const volScalarField& rho = rhoTP();
    const surfaceScalarField alphaRhoPhi = this->alphaRhoPhi_*fvc::interpolate(rho);
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    #if OFVERSION >= 900
        const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
        const Foam::fvConstraints& fvConstraints
        (
            Foam::fvConstraints::New(this->mesh_)
        );
    #else 
        fv::options& fvOptions(fv::options::New(this->mesh_));
    #endif

    BasicMomentumTransportModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
            #if OFVERSION >= 900
          + fvModels.source(alpha, rho, omega_)
            #else 
          + fvOptions(alpha, rho, omega_)
            #endif
        );

        omegaEqn.ref().relax();
        #if OFVERSION >= 900
            fvConstraints.constrain(omegaEqn.ref());
            omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
            solve(omegaEqn);
            fvConstraints.constrain(omega_);
        #else 
            fvOptions.constrain(omegaEqn.ref());
            omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
            solve(omegaEqn);
            fvOptions.correct(omega_);
        #endif
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_)
      + kSource()
        #if OFVERSION >= 900
      + fvModels.source(alpha, rho, k_)
        #else 
      + fvOptions(alpha, rho, k_)
        #endif
    );

    kEqn.ref().relax();
    #if OFVERSION >= 900
        fvConstraints.constrain(kEqn.ref());
        solve(kEqn);
        fvConstraints.constrain(k_);
    #else 
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(k_);
    #endif
    bound(k_, this->kMin_);

    correctNut(S2, F23);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
