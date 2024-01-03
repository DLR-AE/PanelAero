# Imports from python
import pickle
# Imports from loadskernel
from loadskernel import build_aero_functions
# Import from this repository
from panelaero import VLM, DLM


class NewModel():

    def build_aerogrid(self, filename):
        self.aerogrid = build_aero_functions.build_aerogrid(filename, method_caero='CAERO1')


# Build the aerogrid
model = NewModel()
model.build_aerogrid('../../example/simplewing.CAERO1')

with open('./simplewing_aerogrid.pickle', 'wb') as fid:
    pickle.dump(model.aerogrid, fid)

# VLM
Qjj = VLM.calc_Qjj(model.aerogrid, Ma=0.0)
with open('./simplewing_VLM_Ma00.pickle', 'wb') as fid:
    pickle.dump(Qjj, fid)

Qjj = VLM.calc_Qjj(model.aerogrid, Ma=0.3)
with open('./simplewing_VLM_Ma03.pickle', 'wb') as fid:
    pickle.dump(Qjj, fid)

# DLM parabolic
Qjj = DLM.calc_Qjj(model.aerogrid, Ma=0.0, k=0.2)
with open('./simplewing_DLM_Ma00_k02.pickle', 'wb') as fid:
    pickle.dump(Qjj, fid)

Qjj = DLM.calc_Qjj(model.aerogrid, Ma=0.3, k=0.2)
with open('./simplewing_DLM_Ma03_k02.pickle', 'wb') as fid:
    pickle.dump(Qjj, fid)


# DLM quartic
Qjj = DLM.calc_Qjj(model.aerogrid, Ma=0.0, k=0.2, method='quartic')
with open('./simplewing_DLM_Ma00_k02_quartic.pickle', 'wb') as fid:
    pickle.dump(Qjj, fid)

Qjj = DLM.calc_Qjj(model.aerogrid, Ma=0.3, k=0.2, method='quartic')
with open('./simplewing_DLM_Ma00_k02_quartic.pickle', 'wb') as fid:
    pickle.dump(Qjj, fid)

print('Done.')
