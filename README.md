# An Implementation of the Vortex Lattice and the Doublet Lattice Method
A Vortex Lattice Method (VLM) and a Doublet Lattice Method (DLM) is implemented in Python. The aerodynamic influence matrices (AICs) obtained from this implementation are validated with respect to MSC.Nastran for both the parabolic and the quartic integration schemes of the DLM. The test cases include dihedral and sweep of the main wing, wing-empennage configurations with the horizontal tail planar to the main wing, co-planar and further away (e.g. T-tail) and with/without a vertical tail. The test cases have been inspected at different mach numbers and reduced frequencies. For all tested aircraft configurations the results were found to be equivalent to MSC.Nastran in a numerical sense. Using panels which are misaligned in y-direction provokes differences and errors.

# Reference
Voß, A., “An Implementation of the Vortex Lattice and the Doublet Lattice Method,” Institut für Aeroelastik, Deutsches Zentrum für Luft- und Raumfahrt, Göttingen, Germany, Technical Report DLR-IB-AE-GO-2020-137, Oktober 2020, https://elib.dlr.de/136536/.

If you use this software for your scientific work, we kindly ask you to include a reference to this report in your publications. Thank you!

# Installation & Use
## User installation 
To install everything as a python package, including dependencies:

```
pip install --user git+https://github.com/DLR-AE/PanelAero.git
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
git clone https://github.com/DLR-AE/PanelAero.git
cd ./panel-aero
pip install --user -e . 
```

# License
This software is developed for scientific applications and is delivered as open source without any liability (BSD 3-Clause, please see the [license](LICENSE) for details). For every new aircraft, a validation against test data and/or other simulation tools is highly recommended and in the responsibility of the user. 

If you use this software for your scientific work, we kindly ask you to include a reference to the technical report (see above) in your publications. Thank you!

# Feedback & Support
Note that this is a scientific software for users with a background in aerospace engineering and with a good understanding and experience in aeroelasticity. If you know what you are doing - go ahead and have fun! If you need specific help or assistence, we offer commerical support:
- Development of additional, proprietary features
- Consulting & Training courses
- Service & Support

We are interested in partnerships from both industry and academia, so feel free to contact us (arne.voss@dlr.de).

If you discoverd an obvious bug, please open an [issue](https://github.com/DLR-AE/PanelAero/issues). In case you already know how to fix it, please provide your feedback via merge requests. For details, please see the [instructions](CONTRIBUTING.md) on how to provide a contribution or contact arne.voss@dlr.de if you need any assistance with that.

# Continuous Integration
Status of the (internal) DLR GitLab continuous integration pipelines:

Master branch [![pipeline status](https://gitlab.dlr.de/loads-kernel/panel-aero/badges/master/pipeline.svg)](https://gitlab.dlr.de/loads-kernel/panel-aero/-/commits/master)

Development branch [![pipeline status](https://gitlab.dlr.de/loads-kernel/panel-aero/badges/devel/pipeline.svg)](https://gitlab.dlr.de/loads-kernel/panel-aero/-/commits/devel)

Test coverage [![coverage](https://gitlab.dlr.de/loads-kernel/panel-aero/badges/master/coverage.svg)](https://loads-kernel.pages.gitlab.dlr.de/panel-aero/coverage/)