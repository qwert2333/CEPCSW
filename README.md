# [CEPCSW](https://cepc.github.io/CEPCSW/)

[![Build Status](https://www.travis-ci.com/cepc/CEPCSW.svg?branch=master)](https://www.travis-ci.com/cepc/CEPCSW)
[![CI](https://github.com/cepc/CEPCSW/workflows/CI/badge.svg?branch=master)](https://github.com/cepc/CEPCSW/actions)

CEPC offline software prototype based on [Key4hep](https://github.com/key4hep).

## Quick start

SSH to lxslc7 (CentOS 7).

Before run following commands, please make sure you setup the CVMFS:

```
$ git clone git@github.com:cepc/CEPCSW.git
$ cd CEPCSW
$ git checkout dev-EcalRec-v05 # branch name
$ source setup.sh
$ ./build.sh
$ ./run.sh Examples/options/helloalg.py
```

## Packages

* Examples: For new comers and users

* Detector: Geometry

* Generator: Physics Generator

* Simulation: Detector Simulation

* Digitization: Digitization

* Reconstruction: Reconstruction


## Conventions for collections
Keep the collection names compatible between the prototype and the existing CEPC software.

* MCParticle
* VXDCollection
* SITCollection
* TPCCollection
* SETCollection


## Update in branch dev-EcalRec-v05: CRD Ecal Reconstruction
  * Focus on hadron reconstruction in Ecal.
  * NOTE: This branch has been abandoned. Arbor connection in 1D is accomplished, but matching them in 2D is hard to do. Please go to next version.  

