# Continuous Phase Area Detection Solver/PostProcess library for OpenFOAM

This is a very small library for the detection of continuos phase areas (cpas) in multiphase 
VoF Simulations in OpenFOAM (OF). It may be used for Droplet, Bubble or Film Tracking 
in transient VoF Simulations in OF.
It can be used during runtime as addition to an existing two phase VoF solver, or as function in a custom post processing utility.

It was developed from me for the automated detection and localization of simulated droplets for the investigation
of the dripping behavior during the electroslag remelting (ESR) process.

The following properties are written to the report file(s):
  - x,y,z location
  - mass
  - volume 
  - alpha_max
  - alpha_mean (volume averaged)

You may choose between detection between cells faces:

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

reportFileName  "cpa.txt"; // filename suffix
targetPhase     1;         // alpha1
alphaMinThreshold 0.2;     // min. threshold for continuous area detection
detectOverCellNodes on;    // continuous area detection over cells-edges(nodes)=on or cells-faces=off
```

### Parallel considerations
This library also works for parallel cases, but there is no consolidation of cpas which belong together, but are located on different processors, during runtime. Instead each processor writes it's own results file, where necessary parallel information is given. This was done for ease of implementation and to reduce communication overhead during runtime. The different files can be reconstructed using the Julia scripts, which can be found under ./eval.

## Evaluation of results
Under ./eval there a some exemplary julia script and functions which should make the evaluation of results quite simple.