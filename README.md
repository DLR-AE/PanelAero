# An Implementation of the Vortex Lattice and the Doublet Lattice Method
A Vortex Lattice Method (VLM) and a Doublet Lattice Method (DLM) is implemented in Python. The aerodynamic influence matrices (AICs) obtained from this implementation are validated with respect to MSC.Nastran for both the parabolic and the quartic integration schemes of the DLM. The test cases include dihedral and sweep of the main wing, wing-empennage configurations with the horizontal tail planar to the main wing, co-planar and further away (e.g. T-tail) and with/without a vertical tail. The test cases have been inspected at different mach numbers and reduced frequencies. For all tested aircraft configurations the results were found to be equivalent to MSC.Nastran in a numerical sense. Using panels which are misaligned in y-direction provokes differences and errors.

# Reference
Voß, A., “An Implementation of the Vortex Lattice and the Doublet Lattice Method,” Institut für Aeroelastik, Deutsches Zentrum für Luft- und Raumfahrt, Göttingen, Germany, Technical Report DLR-IB-AE-GO-2020-137, Oktober 2020, https://elib.dlr.de/136536/.

# Continuous Integration

Master branch [![pipeline status](https://gitlab.dlr.de/loads-kernel/panel-aero/badges/master/pipeline.svg)](https://gitlab.dlr.de/loads-kernel/panel-aero/-/commits/master)

Development branch [![pipeline status](https://gitlab.dlr.de/loads-kernel/panel-aero/badges/devel/pipeline.svg)](https://gitlab.dlr.de/loads-kernel/panel-aero/-/commits/devel)

Test coverage [![coverage](https://gitlab.dlr.de/loads-kernel/panel-aero/badges/master/coverage.svg)](https://loads-kernel.pages.gitlab.dlr.de/panel-aero/coverage/)

# Installation & Use
## User installation 
To install everything as a python package, including dependencies:

```
pip install --user git+https://gitlab.dlr.de/loads-kernel/panel-aero.git 
```
## How can I use it?

In Python, you can import the VLM or the DLM as shown below. For further details, please see the example section.

```
import panelaero.VLM as VLM
import panelaero.DLM as DLM
```

## Developer installation 
As above, but with access to the code (keep the code where it is so that you can explore and modify):

```
git clone https://gitlab.dlr.de/loads-kernel/panel-aero.git 
cd ./panel-aero
pip install --user -e . 
```



## Feedback
Please provide your feedback via merge requests (please see [CONTRIBUTING.md](CONTRIBUTING.md)
for details) or contact Arne Voß, arne.voss@dlr.de.
