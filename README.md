# Continuous Phase Area Detection Solver/PostProcess library for OpenFOAM

This is a very small library for the detection of continuos phase areas (cpas) in multiphase 
VoF Simulations in OpenFOAM (OF). It may be used for droplet, bubble or film Tracking 
in transient VoF Simulations in OF.
It can be used during runtime as addition to an existing two phase VoF solver, or as function in a custom post processing utility.

It was developed from me for the automated detection and localization of simulated droplets for the investigation
of the dripping behavior during the electroslag remelting (ESR) process.

The following properties are written to the report file(s) for each continous phase area
for each timestep:

  - no of cells
  - x,y,z location
  - mass
  - volume 
  - alpha_max
  - alpha_mean (volume averaged)
  - connected boundaries

You may choose between detection between cells faces, cell edges and dell points, i.e.:

![GitHub Logo](/images/dropletDetection_cn_vs_cf.png)


## How to use
The library was written and tested for OpenFOAM 1906+. I is reasonable fast, but I would neither call it highly sophisticated nor very well tested, so use it on your own risk. It can be used as addition to a solver or as a post processing utility on existing results (see ./examples/postprocessing).

### Solver Modification
The continuous phase area detection should work for all **two** phase VoF solvers, where the `volScalarField`(s) `alpha1` and `alpha2` are available.

First, copy cpa.C, cpa.H, cpad.C, cpad.H, getLocalToGlobalFaceArray.H to the solvers folder
Then, modify the solvers C File to include the following changes:

```c++
#include <cpad.H>
```

Append the following code at the end of the runTime loop:

```c++
while (runTime.run())
{
    .
    .
    .


    #include "getLocalToGlobalFaceArray.h"

    FoamCpad::cpad cpas(mixture, runTime);
    cpas.appendTimeStepReport(
                                runTime.value(),
                                procFaceToGlobalFaceList
                             );
}
```

Theoretically, #include "getLocalToGlobalFaceArray.h" could be moved out of runTime loop, for unchanging meshes and getting little better performance, but in my experience the cpad detection time was negligible compared to calculation time anyhow, but that may be not the case for different cases and meshes...



At least, change the Makefiles (files and options) of the solver to include the following modifications. Assuming the cpad source files are located in the parent directory of the solver you are going to modify:

```
/... 
    | - cpad/
            | Make/
                  | files
                  | options
            | Cpa/
            | cpad.C
            | cpad.H
    | - Solver/
            | Make/
                  | files
                  | options
            | solvercpad.c
            | ...
```
Makefile modifications:

files:  
*(You probably want to rename the solver to keep the standard unmodified solver...)*
```Makefile
solvercpad.C

EXE = $(FOAM_USER_APPBIN)/solvercpad
```

options:

At least the following modifications have to be done. Due to latest changes to support solver, which either use twoPhaseMixutreThermo, twoPhaseMixture or immiscibleIncompressibleTwoPhaseMixture solver you have one of the following options:

1. You have to make sure that your modified solver Make/options file includes all libraries from the cpad library. 
2. Remove the unecessary mixtureTypes from cpad.H (includes and constructors) and recompile the library.

But **at least** the following changes will be necessary:
```Makefile
EXE_INC = \
    .
    .
    .
    -I../cpad \
    .
    .
    .
    -I$(LIB_SRC)/transportModels/cpad/lnInclude

EXE_LIBS = \
    .
    .
    .
    -L$(FOAM_USER_LIBBIN) \
    -lcpad \
```

Compile and test!


### The system/cpad Dictionary
For settings of either post process utility or modified solver there it is expected to have a config file named "cpad" under the cases system/ directory. 

```c++
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.0                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      cpad;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

reportFileName      "cpa.txt";  // filename suffix
targetPhase         1;          // alpha1
alphaMinThreshold   0.2;        // min. threshold for continuous area detection
detectOver          "faces";    // continuous area detection over: "faces" "edges" or 
                                // "points"
```

### Parallel considerations
This library also works for parallel cases, but there is no consolidation of cpas which belong together, but are located on different processors, during runtime. Instead each processor writes it's own results file, where necessary parallel information is given. 
This was done for ease of implementation and to reduce communication overhead during runtime. The different files can be reconstructed using the Julia scripts, which can be found under ./eval.

**Currently any other detections other than "over faces", is not working across processor boundary faces.** So if ``detectOver`` "edges" or "points" is choosen, there may, under certain circumstances, be deviations between serial and parallel cpa evaluations. 


## Evaluation of results
Under ./eval there a some exemplary julia script and functions which should make the evaluation of results quite simple.


## Disclaimer

*The project in which the work on this solver was fundet is over by now, so I wont be able to make many enhancements in the near future myself. We hope there will be future projects to make more enhancements to these solver collection. Nevertheless if there may be any pull request I will try to have a look ...*

## Acknowledgment
The authors gratefully acknowledge the support of the
German Research Foundation (Deutsche Forschungsgemeinschaft – DFG) during the DFG project PF 394/24-1 (2016-2019).