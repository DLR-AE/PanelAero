"""
This is a short demonstration on how to use the DLM at the example of a simple wing.
Note that functionalities of the Loads Kernel are used for
- building the aerogrid from a CAERO card and for
- plotting the output.
Presumably, you will replace these two steps during the integration of the DLM into your own software anyway.
"""
# Imports from python
import numpy as np
from matplotlib import pyplot as plt
import sys
# Import from this repository
import panelaero.DLM as DLM
# Imports from loadskernel
# Here you add the location of the Loads Kernel
sys.path.append("../../loads-kernel")
from loadskernel import build_aero_functions, plotting_extra

# Geometrie
aerogrid = build_aero_functions.build_aerogrid('./simplewing.CAERO1', method_caero = 'CAERO1')
# DLM
Qjj = DLM.calc_Qjj(aerogrid, Ma=0.0, k=0.1) # with k = omega/U, the "classical" definition, not Nastran definition!

# Anströmung
Vtas = 25.0
q_dyn = 1.225/2.0 * Vtas**2
wj = np.ones(aerogrid['n']) * 1.0/Vtas # downwash von 1.0 m/s auf jedes Panel, skaliert mit Vtas
# Beispielrechnung
cp = Qjj.dot(wj)
Fxyz = q_dyn * aerogrid['N'].T * aerogrid['A'] * Qjj.dot(wj)

# Plot von Real und Imaginärteil
ax1 = plotting_extra.DetailedPlots.plot_aerogrid(plotting_extra.DetailedPlots, aerogrid, cp.real, colormap = 'plasma')
ax2 = plotting_extra.DetailedPlots.plot_aerogrid(plotting_extra.DetailedPlots, aerogrid, cp.imag, colormap = 'plasma')
plt.show()

print('Done.')