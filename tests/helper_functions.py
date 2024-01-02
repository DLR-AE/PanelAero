import numpy as np


class HelperFunctions(object):

    def compare_AICs(self, AIC_a, AIC_b, n):
        result_allclose = np.allclose(AIC_a, AIC_b)
        print('Numerically equal? {}'.format(result_allclose))
        if not result_allclose:
            # How large is the difference?
            print('Sum of differences = {}'.format(str(np.sum(AIC_a - AIC_b))))
            # Are the difference in the real or in the imaginary part?
            M_real = []
            M_imag = []
            for i_row in range(n):
                M_real.append(np.dot(AIC_b[i_row, :].real, AIC_a[i_row, :].real) / np.linalg.norm(AIC_b[i_row, :].real)
                              / np.linalg.norm(AIC_a[i_row, :].real))
                M_imag.append(np.dot(AIC_b[i_row, :].imag, AIC_a[i_row, :].imag) / np.linalg.norm(AIC_b[i_row, :].imag)
                              / np.linalg.norm(AIC_a[i_row, :].imag))
            print('m_real = {}, m_imag = {}'.format(np.mean(M_real), np.mean(M_imag)))
        return result_allclose

#             # Some plots to figure out the location/source of the differences.
#             plot = DetailedPlots(self.jcl, self)
#             ax = plot.plot_aerogrid(self.aerogrid, M_real, 'bwr', -1.0, 1.0)
#             ax = plot.plot_aerogrid(self.aerogrid, M_imag, 'bwr', -1.0, 1.0)
#
#             wj = np.ones(self.aerogrid['n'])*5.0/180.0*np.pi
#             ax = plot.plot_aerogrid(self.aerogrid, AIC_a.dot(wj).real - AIC_b.dot(wj).real, 'bwr', -0.1, 0.1)
#             ax = plot.plot_aerogrid(self.aerogrid, AIC_a.dot(wj).imag - AIC_b.dot(wj).imag, 'bwr', -0.1, 0.1)
#             ax = plot.plot_aerogrid(self.aerogrid, AIC_a.dot(wj).real, 'plasma')
#             ax = plot.plot_aerogrid(self.aerogrid, AIC_b.dot(wj).real, 'plasma')
#             plt.show()
