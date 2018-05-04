
from __future__ import division, print_function
from itertools import product
import numpy as np

from .utils import safe_inverse


class Mixer(object):
    """DAgostini Bayesian Mixer Class
    """
    def __init__(self, error_type='ACM', MCTables=None, data=None,
                 data_err=None, efficiencies=None, efficiencies_err=None,
                 response=None, response_err=None):
        # Input validation
        if len(data) != response.shape[0]:
            err_msg = ('Inconsistent number of effect bins. Observed data '
                       'has {} effect bins, while response matrix has {} '
                       'effect bins.'.format(len(data), response.shape[0]))
            raise ValueError(err_msg)

        # Normalized P(E|C)
        self.pec = response
        self.cEff = efficiencies
        self.NEobs = data
        self.error_type = error_type

        # Number of Cause and Effect Bins
        dims = self.pec.shape
        self.cbins = dims[1]
        self.ebins = dims[0]
        self.cEff_inv = safe_inverse(self.cEff)
        # Mixing Matrix
        self.Mij = np.zeros(dims)

        self.cov = CovarianceMatrix(error_type=error_type,
                                    data=data,
                                    data_err=data_err,
                                    efficiencies=efficiencies,
                                    efficiencies_err=efficiencies_err,
                                    response=response,
                                    response_err=response_err)

    def get_cov(self):
        """Covariance Matrix
        """
        cvm = self.cov.get_cov()
        return cvm

    def get_stat_err(self):
        """Statistical Errors
        """
        cvm = self.cov.getVc0()
        err = np.sqrt(cvm.diagonal())
        return err

    def get_MC_err(self):
        """MC (Systematic) Errors
        """
        cvm = self.cov.getVc1()
        err = np.sqrt(cvm.diagonal())
        return err

    def smear(self, n_c):
        """Smear Calculation - the Bayesian Mixer!

        Only needs input prior, n_c, to unfold
        """
        if len(n_c) != self.cbins:
            err_msg = ('Trying to unfold with the wrong number of causes. '
                       'Response matrix has {} cause bins, while prior '
                       'has {} cause bins.'.format(self.cbins, len(n_c)))
            raise ValueError(err_msg)

        ebins = self.ebins
        cbins = self.cbins

        # Bayesian Normalization Term (denominator)
        f_norm = np.dot(self.pec, n_c)
        f_inv = safe_inverse(f_norm)

        # Unfolding (Mij) Matrix at current step
        Mij = np.zeros(self.Mij.shape)

        n_c_eff = n_c * self.cEff_inv
        for ti, ej in product(range(0, cbins), range(0, ebins)):
            Mij[ej, ti] = self.pec[ej, ti] * n_c_eff[ti] * f_inv[ej]

        # Estimate cause distribution via Mij
        n_c_update = np.dot(self.NEobs, Mij)

        # The status quo
        self.Mij = Mij.copy()
        self.cov.set_current_state(self.Mij, f_norm, n_c_update, n_c)

        return n_c_update


class CovarianceMatrix(object):
    """Covariance Base Class

    All Covariance Matrix Code either transcribed from Adye's RooUnfoldBayes.cxx
    (http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html#dev),
    or based on derivations presented in Unfolding reference section 'Unfolding Uncertainties'
    """
    def __init__(self, error_type='ACM', MCTables=None, data=None, data_err=None,
                 efficiencies=None, efficiencies_err=None, response=None,
                 response_err=None):

        self.error_type = error_type
        # Normalized P(E|C)
        self.pec = response
        self.pec_err = response_err
        # Errors on P(E|C), for Vc1
        self.pec_cov_type = "Multinomial"
        # self.pec_cov_type = "Poisson"
        self.PecCov = self.SetPecCovType(self.pec_cov_type)
        self.cEff = efficiencies
        self.cEff_inv = safe_inverse(self.cEff)
        # Effective number of sim events
        efficiencies_err_inv = safe_inverse(efficiencies_err)
        NCmc = (efficiencies * efficiencies_err_inv)**2
        self.NCmc = NCmc
        self.NEobs = data
        self.NEobs_err = data_err
        # Total number of observed effects
        self.nobs = np.sum(self.NEobs)

        # Number of cause and effect eins
        dims = self.pec.shape
        self.cbins = dims[1]
        self.ebins = dims[0]
        # Mixing matrix
        self.Mij = np.zeros(dims)
        # Flag to propagate derivative
        self.ErrorPropFlag = self.SetErrorPropFlag(self.error_type)
        # Adye propagating derivative
        self.dcdn = np.zeros(dims)
        self.dcdP = np.zeros((self.cbins, self.cbins * self.ebins))
        # Counter for number of iterations
        self.counter = 0

    def set_current_state(self, Mij, f_norm, n_c, n_c_prev):
        """Set the Current State of dcdn and dcdP
        """
        # For ease of typing
        ebins = self.ebins
        cbins = self.cbins

        # First set Mij
        self.Mij = Mij.copy()

        # D'Agostini Form (and/or First Term of Adye)
        dcdn = Mij.copy()
        dcdP = np.zeros(self.dcdP.shape)
        f_inv = safe_inverse(f_norm)

        NE_F_R = self.NEobs * f_inv
        # (ti, ec_j + tk) elements
        for ej, ti, tk in product(range(0, ebins), range(0, cbins), range(0, cbins)):
            ec_j = ej * cbins
            dcdP[ti, ec_j+tk] = -NE_F_R[ej] * Mij[ej, ti] * n_c_prev[tk]
        # (ti, ec_j + ti) elements
        for ej, ti in product(range(0, ebins), range(0, cbins)):
            ec_j = ej * cbins
            dcdP[ti, ec_j+ti] += (n_c_prev[ti] * NE_F_R[ej] - n_c[ti]) * self.cEff_inv[ti]

        # Adye propagation corrections
        if self.ErrorPropFlag and self.counter > 0:
            # Get previous derivatives
            dcdn_prev = self.dcdn.copy()
            dcdP_prev = self.dcdP.copy()

            n_c_prev_inv = safe_inverse(n_c_prev)
            # Ratio of updated n_c to n_c_prev
            nc_r = n_c * n_c_prev_inv
            # Efficiency ratio of n_c_prev
            e_r = self.cEff * n_c_prev_inv

            # Calculate extra dcdn terms
            M1 = dcdn_prev * nc_r
            M2 = -Mij * e_r
            M2 = M2.T * self.NEobs
            M3 = np.dot(M2, dcdn_prev)
            dcdn += np.dot(Mij, M3)
            dcdn += M1

            # Calculate extra dcdP terms
            # My version (from unfolding doc)
            At = Mij.T * self.NEobs
            B = Mij * e_r
            C = np.dot(At, B)
            dcdP_Upd = np.dot(C, dcdP_prev)
            dcdP += (dcdP_prev.T * nc_r).T - dcdP_Upd

        # Set current derivative matrices
        self.dcdn = dcdn.copy()
        self.dcdP = dcdP.copy()
        # On to next iteration
        self.counter += 1

    def getVcd(self):
        """Get Covariance Matrix of N(E), ie from Observed Effects
        """
        Vcd = np.diag(self.NEobs_err**2)
        return Vcd

    def getVc0(self):
        """Get full Vc0 (data) contribution to cov matrix
        """
        # Get derivative
        dcdn = self.dcdn.copy()
        # Get NObs covariance
        Vcd = self.getVcd()
        # Set data covariance
        Vc0 = dcdn.T.dot(Vcd).dot(dcdn)

        return Vc0

    def getVcPP(self):
        """Get Covariance Matrix of P(E|C), ie from MC
        """
        cbins = self.cbins
        ebins = self.ebins

        CovPP = np.zeros((cbins * ebins, cbins * ebins))

        # Poisson covariance matrix
        if self.PecCov == 1:
            for ej, ti in product(range(0, ebins), range(0, cbins)):
                ejc = ej * cbins
                CovPP[ejc+ti, ejc+ti] = self.pec_err[ej, ti]**2
        # Multinomial covariance matrix
        elif self.PecCov == 2:
            NC_inv = safe_inverse(self.NCmc)
            for ej, ti in product(range(0, ebins), range(0, cbins)):
                ejc = ej * cbins
                # Don't go looping if zeros are present :)
                if (self.pec[ej, ti] > 0 and NC_inv[ti] > 0):
                    CovPP[ejc+ti, ejc+ti] = NC_inv[ti] * self.pec[ej, ti] * (1 - self.pec[ej, ti])
                    cov1 = -NC_inv[ti] * self.pec[ej, ti]
                    for ek in range(ej+1, ebins):
                        ekc = ek * cbins
                        if ejc+ti == ekc+ti or ek == ej:
                            continue
                        cov = cov1 * self.pec[ek, ti]
                        CovPP[ejc+ti, ekc+ti] = cov
                        CovPP[ekc+ti, ejc+ti] = cov

        return CovPP

    def getVc1(self):
        """Get full Vc1 (MC) contribution to cov matrix
        """
        # Get NObs covariance
        CovPP = self.getVcPP()
        # Get derivative
        dcdP = self.dcdP
        # Set MC covariance
        Vc1 = dcdP.dot(CovPP).dot(dcdP.T)
        return Vc1

    def get_cov(self):
        """Get full covariance matrix
        """
        # Get Contributions from Vc0 and Vc1
        Vc0 = self.getVc0()
        Vc1 = self.getVc1()
        # Full Covariance Matrix
        Vc = Vc0 + Vc1
        return Vc

    def SetErrorPropFlag(self, func):
        """Key to choose whether to propagate errors at each iteration
        """
        # ACM propagate errors as Adye
        # DCM propagate errors as DAgostini
        ErrorPropOptionsKey = {
            "ACM": True,
            "DCM": False,
            }
        if func in ErrorPropOptionsKey:
            return ErrorPropOptionsKey[func]
        else:
            err_msg = ('Invalid error propagation type entered, {}, must '
                       'be in {}'.format(func, ErrorPropOptionsKey.keys()))
            raise ValueError(err_msg)

    def SetPecCovType(self, func):
        """Key for Selecting PEC Covariance Matrix Type
        """
        PecCovOptionsKey = {
            "Poisson": 1,
            "Multinomial": 2
            }
        if func in PecCovOptionsKey:
            return PecCovOptionsKey[func]
        else:
            err_msg = ('Invalid covariance calculation type entered, {}, must '
                       'be in {}'.format(func, PecCovOptionsKey.keys()))
            raise ValueError(err_msg)
