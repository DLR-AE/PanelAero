"""
This is a short demonstration on how to use the DLM at the example of a simple wing.
"""
# Imports from python
import numpy as np

# Import from this repository
from panelaero import DLM
from example.build_aeromodel import AeroModel
from example.plotting import DetailedPlots

# build a model that includes the aerodynmic grid
model = AeroModel('./simplewing.CAERO1')
model.build_aerogrid()

# run the DLM
# with k = omega/U, the "classical" definition, not Nastran definition!
Qjj = DLM.calc_Qjj(model.aerogrid, Ma=0.0, k=0.1)

# set-up some generic onflow
Vtas = 25.0
q_dyn = 1.225 / 2.0 * Vtas ** 2
# downwash of 1.0 m/s on every panel, scaled with Vtas
wj = np.ones(model.aerogrid['n']) * 1.0 / Vtas
# calculate the pressure coefficent distribution
cp = Qjj.dot(wj)

# Plot of real and imaginary parts, note that these plots use mayavi and tvtk
# (possibly not installed by default and to keep the list of dependencies small).
plots = DetailedPlots(model=model)
plots.plot_aerogrid(cp.real)
plots.plot_aerogrid(cp.imag)

# a force vector is calculated like this
Fxyz = q_dyn * model.aerogrid['N'].T * model.aerogrid['A'] * Qjj.dot(wj)

print('Done.')
