import pickle

from panelaero import VLM, DLM
from tests.helper_functions import HelperFunctions


class TestSimplewing(HelperFunctions):
    # Load geometry
    with open('./tests/reference_data/simplewing_aerogrid.pickle', 'rb') as fid:
        aerogrid = pickle.load(fid)

    def test_VLM_Ma00(self):
        # Calculate the AIC matrix
        Qjj = VLM.calc_Qjj(self.aerogrid, Ma=0.0)
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_VLM_Ma00.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjj, reference_data, self.aerogrid['n']), "AIC does NOT match reference"

    def test_VLM_Ma03(self):
        # Calculate the AIC matrix
        Qjj = VLM.calc_Qjj(self.aerogrid, Ma=0.3)
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_VLM_Ma03.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjj, reference_data, self.aerogrid['n']), "AIC does NOT match reference"

    def test_VLM_sequence_of_mach_numbers(self):
        # Calculate the AIC matrix
        Qjjs = VLM.calc_Qjjs(self.aerogrid, Ma=[0.7, 0.3])
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_VLM_Ma03.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjjs[0][1, :, :], reference_data[0], self.aerogrid['n']), "AIC does NOT match reference"

    def test_DLM_Ma00_k02(self):
        # Calculate the AIC matrix
        Qjj = DLM.calc_Qjj(self.aerogrid, Ma=0.0, k=0.2)
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_DLM_Ma00_k02.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjj, reference_data, self.aerogrid['n']), "AIC does NOT match reference"

    def test_DLM_Ma03_k02(self):
        # Calculate the AIC matrix
        Qjj = DLM.calc_Qjj(self.aerogrid, Ma=0.3, k=0.2)
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_DLM_Ma03_k02.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjj, reference_data, self.aerogrid['n']), "AIC does NOT match reference"

    def test_DLM_Ma00_k02_quartic(self):
        # Calculate the AIC matrix
        Qjj = DLM.calc_Qjj(self.aerogrid, Ma=0.0, k=0.2, method='quartic')
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_DLM_Ma00_k02_quartic.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjj, reference_data, self.aerogrid['n']), "AIC does NOT match reference"

    def test_DLM_Ma03_k02_quartic(self):
        # Calculate the AIC matrix
        Qjj = DLM.calc_Qjj(self.aerogrid, Ma=0.3, k=0.2, method='quartic')
        # Do the comparison
        print('Comparing AIC with reference')
        with open('./tests/reference_data/simplewing_DLM_Ma03_k02_quartic.pickle', 'rb') as fid:
            reference_data = pickle.load(fid)
        assert self.compare_AICs(Qjj, reference_data, self.aerogrid['n']), "AIC does NOT match reference"
