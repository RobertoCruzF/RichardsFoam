/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Application
    RichardsFoam

Description
    Transient solver for flow in unsaturated porous media
    With chord slope formulation of the Richards equation.
    van Genuchten laws for unsaturated hydraulic properties parametrisation
    Global computation of the convergence criterium
    Adaptative time stepping with a stabilisation procedure
    NB 1: use backward scheme for time discretisation
    NB 2: use only mesh with constant cell volumes

References
   version 0.0 (develloped with OpenFOAM 2.0.1)
   Details may be found in:
   Orgogozo, L., Renon, N., Soulaine, C., Hénon, F., Tomer, S.K., Labat, D., 
    Pokrovsky, O.S., Sekhar, M., Ababou, R., Quintard, M., Submitted. 
   Mechanistic modelling of water fluxes at the watershed scale: An open source 
    massively parallel solver for Richards equation.
   Submitted to Computer Physics Communications.
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    #include "readTimeControls.H"

    #include "readPicardControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// initialisation of the scalar containing the number of mesh cells. Note the use of gSum instead of sum.
	double nbMesh;
        nbMesh = gSum(usf);

// initialisation of the scalar containing the residual for the exit test of the Picard loop.
        double crit;
        crit=0.;

// initialisation of the token which counts the number of Picard iteraion for the adaptive time step procedure.
        int currentPicard;
	currentPicard = nIterPicard-3;

// initialisation of the token which counts the number of Stabilisation cycle for the stabilisation of the adaptive time step procedure.
        int sc;
	sc = 0;

// initialisation of the field of altitudes.     
    volVectorField positionVector=mesh.C();
    volScalarField z=positionVector.component(vector::Z);

// initialisation of the intermediate fields for the Picard loop.
    volScalarField psi_tmp = psi;
    volScalarField psim1=psi;

// initialisation of the varying transport properties for the Picard loop.
    volScalarField thtil=0.5*((1+sign(psi))+(1-sign(psi))*pow((1+pow(mag(alpha*psi),n)),-(1-(1/n))));
    volScalarField thtil_tmp=0.5*((1+sign(psi_tmp))+(1-sign(psi_tmp))*pow((1+pow(mag(alpha*psi_tmp),n)),-(1-(1/n))));
    volScalarField Krel=0.5*((1+sign(psi))*K+(1-sign(psi))*K*pow(thtil,0.5)*pow((1-pow((1-pow(thtil,(n/(n-1)))),(1-(1/n)))),2)); 
    volScalarField Crel=S+0.5*((1-sign(psi))*((thetas-thetar)*(thtil-thtil_tmp)*(1./((usf*pos(psi-psi_tmp)*pos(psi_tmp-psi))+psi-psi_tmp))));

// initialisation of the gravity term.
    volVectorField gradk=fvc::grad(Krel);
    volScalarField gradkz=gradk.component(vector::Z);

//initialisation of the velocity field.
    U=-Krel*((fvc::grad(psi))+vuz);
         
    Info<< "\nStarting time loop\n" << endl;

// starting of the time loop.
    while (runTime.loop())
    {

// time step control operations.
        #include "readTimeControls.H"
        #include "setDeltaT.H" 

        Info<< "Time = " << runTime.timeName() << nl << endl;

// beginning of the stabilisation loop for the stabilised adaptive time step procedure.
        for (int cyc=0; cyc<nMaxCycle; cyc++)
            {

// beginning of the Picard loop.
        for (int pic=0; pic<nIterPicard; pic++)
             {   
    
              psim1=psi;

	      psi = psi_tmp;

// Resolution of the linear system.
              fvScalarMatrix psiEqn
              (
		          Crel*fvm::ddt(psi)
		          ==
                          fvm::laplacian(Krel,psi,"laplacian(Krel,psi)")
                          +gradkz
              );
			  psiEqn.relax();
			  psiEqn.solve();

// update of the varying transport properties.
              thtil=0.5*((1+sign(psi))+(1-sign(psi))*pow((1+pow(mag(alpha*psi),n)),-(1-(1/n))));
              Krel=0.5*((1+sign(psi))*K+(1-sign(psi))*K*pow(thtil,0.5)*pow((1-pow((1-pow(thtil,(n/(n-1)))),(1-(1/n)))),2));
              Crel=S+0.5*((1-sign(psi))*((thetas-thetar)*(thtil-thtil_tmp)*(1./((usf*pos(psi-psi_tmp)*pos(psi_tmp-psi))+psi-psi_tmp))));      

// update of the gravity term.
              gradk=fvc::grad(Krel);
              gradkz=gradk.component(2);     

// Computation of the field of residuals for the exit test of the Picard loop.            
              err=psi-psim1;

// Computation of the residual for the exit test of the Picard loop. Note the use of gSum instead of sum.
              crit=gSumMag(err)/nbMesh;

// exit test for the Picard loop.         
			 if (crit < precPicard) 
                         {
                         Info<< " Erreur = " << crit 
                             << " Picard = " << pic 
                             << nl << endl;
                         currentPicard=pic;
                         break;
                         }

// end of the Picard loop.
             }

// exit test for the loop associated with the stabilisation cycles for the adaptive time step procedure.
                         if (crit < precPicard) 
                         {
                         break;
                         }
                         else
                         {
                         Info<< "Criterion not reached, restart time loop iteration with a smaller time step / Error = " << crit << nl << endl;
                         runTime.setDeltaT((1/tFact)*runTime.deltaTValue());
                         Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
                         }

// end of the stabilisation cycles loop.              
             }


// warning test in case of convergence failure of the Picard loop.                         
                         if (crit >= precPicard) 
                         {
                         Info<< "Convergence failure/ Error = " << crit << nl << endl;
                         currentPicard=nIterPicard;
                         }

// Final updating of the result fields before going to the next time iteration.
        psi_tmp = psi;
        
        thtil_tmp=0.5*((1+sign(psi_tmp))+(1-sign(psi_tmp))*pow((1+pow(mag(alpha*psi_tmp),n)),-(1-(1/n))));  

        theta=(thetas-thetar)*thtil+thetar;

        U=-Krel*((fvc::grad(psi))+vuz);

// Writting  of the result.
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

// end of the time loop.
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
