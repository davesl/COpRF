# COpRF
Convex Optimized Population Receptive Field (CO-pRF) toolbox

The COpRF toolbox provides an efficient and robust way to estimate population receptive field (pRF) measures from neuroimaging data. The toolbox is an implementation of the linear framework for Convex Optimized Population Receptive Field (COpRF) Mapping described here:

> **Convex Optimized Population Receptive Field (CO-pRF) Mapping.**  
> Slater D., Melie-Garcia L., Adaszewski S., Lutti A., Kay K., Draganski B., Kherif F., 2016.  
> (In preparation)


### *** Release 1.1 - 22/11/2016 ***
This is the current version of the CO-pRF framework and is implemented as a matlab toolbox.

The current version provides:
* Routines to efficiently estimate visual pRF measures from fMRI data
* Option to estimate pRF measures using a choice of two different pRF models - the original pRF model of Dumoulin and Wandell (2008) or the extended CSS-pRF model from Kay et al. (2013).

Future releases will:
* Add additional receptive field shapes available for visual pRF estimations
* Extend the framework to other sensory modalities (e.g. auditory)
* Provde a user interface
* Allow for a user to specify any arbitrary pRF model formula for use across all sensory and imaging modalities

# Installation 

## COpRF Setup

Download the COpRF files and add the download folder and subfolders to the matlab path.

To see the basic setup and running of the toolbox on synthetic data you can try the example script found at:

> COpRF\Examples\COpRF_example1_SyntheticData.m

## External Toolboxes
All required functions to run COpRF are included within the COpRF toolbox package. Here we make reference to the external packed used within the release.

### knkutils
We make use of the MATLAB utilities developed by Kendrick Kay in his knkutils matlab package. The included functions provide complement the builtin matlab fuctions. The package also provides tools for designing and running visual pRF fMRI experiments.

> Code and utilities: https://github.com/kendrickkay/knkutils  
> Additional information: http://cvnlab.net/analyzePRF/

### SPArse Modeling Software (SPAMS)

We make use of the highly efficient optimization routines implemented as part of the SPAMS software package. Complied `mexLasso` functions are included as part of the COpRF toolbox. If you enounter any errors related to the `mexLasso` function when running the toolbox you can download the matlab source code and compile on your system. The SPAMS toolbox can be found here:

> http://spams-devel.gforge.inria.fr/downloads.html
 
