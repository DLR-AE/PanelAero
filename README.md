# An Implementation of the Vortex Lattice and the Doublet Lattice Method
A Vortex Lattice Method (VLM) and a Doublet Lattice Method (DLM) is implemented in Python. The aerodynamic influence matrices (AICs) obtained from this implementation are validated with respect to MSC.Nastran for both the parabolic and the quartic integration schemes of the DLM. The test cases include dihedral and sweep of the main wing, wing-empennage configurations with the horizontal tail planar to the main wing, co-planar and further away (e.g. T-tail) and with/without a vertical tail. The test cases have been inspected at different mach numbers and reduced frequencies. For all tested aircraft configurations the results were found to be equivalent to MSC.Nastran in a numerical sense. Using panels which are misaligned in y-direction provokes differences and errors.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8343766.svg)](https://doi.org/10.5281/zenodo.8343766)

# Theoretical Background & Reference
Voß, A., “An Implementation of the Vortex Lattice and the Doublet Lattice Method,” Institut für Aeroelastik, Deutsches Zentrum für Luft- und Raumfahrt, Göttingen, Germany, Technical Report DLR-IB-AE-GO-2020-137, Oktober 2020, https://elib.dlr.de/136536/.

If you use this software for your scientific work, we kindly ask you to include a reference to this report in your publications. Thank you!

# Installation & Use
## Basic Installation 
Install Panel Aero as a python package with its core dependencies using a package manager (PyPI or Conda):

```
pip install PanelAero
```

or

```
conda install -c conda-forge PanelAero
```

## How can I use it?
In Python, you can import the VLM or the DLM as shown below. For further details, please see the Tutorials section. This is no stand-alone aerodynamic software but is intended to be integrated in other software, for example for loads and aeroelastic analyses.

```
from panelaero import VLM, DLM
```

## Advanced Installation 
As above, but with access to the code (download and keep the code where it is so that you can explore and modify):

```
git clone https://github.com/DLR-AE/PanelAero.git
cd <local_repo_path>
pip install -e . 
```

## Tutorials & Examples
There is a growing number of tutorials based on Jupyter notebooks. You can either have a look at the static html tutorials or use the Jupyter notebooks interactively. For the latter, start a jupyter notebook server, which will open a dashboard in your web browser. Then open one of the *.ipynb notebooks from ./doc/tutorials and walk through the tutorials step-by-step.

[View html tutorials](https://dlr-ae.github.io/PanelAero/tutorials/)

or

```
jupyter notebook
```
Any missing dependencies (probably jupyter and mayavi) can be installed with:

```
pip install -e .[tutorials]
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
Status of the continuous integration pipelines / workflows:

[View test coverage](https://dlr-ae.github.io/PanelAero/coverage/)

Master branch 

[![Regression Tests](https://github.com/DLR-AE/PanelAero/actions/workflows/regression-tests.yml/badge.svg?branch=master)](https://github.com/DLR-AE/PanelAero/actions/workflows/regression-tests.yml)
[![Coding style](https://github.com/DLR-AE/PanelAero/actions/workflows/coding-style.yml/badge.svg?branch=master)](https://github.com/DLR-AE/PanelAero/actions/workflows/coding-style.yml)

Development branch 

[![Regression Tests](https://github.com/DLR-AE/PanelAero/actions/workflows/regression-tests.yml/badge.svg?branch=devel)](https://github.com/DLR-AE/PanelAero/actions/workflows/regression-tests.yml)
[![Coding style](https://github.com/DLR-AE/PanelAero/actions/workflows/coding-style.yml/badge.svg?branch=devel)](https://github.com/DLR-AE/PanelAero/actions/workflows/coding-style.yml)

