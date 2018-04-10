"""
   Covariance Matrix Calculator for D'Agostini
   unfolding method, with modifications to
   error propagation based on method request
   from user input.
"""
from itertools import product
import numpy as np
from .utils import safe_inverse


class CovarianceMatrix(object):
    """Covariance Base Class

    All Covariance Matrix Code either transcribed from Adye's RooUnfoldBayes.cxx
    (http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html#dev),
    or based on derivations presented in Unfolding reference section 'Unfolding Uncertainties'
    """
    def __init__(self, name, MCTables=None, data=None, data_err=None):

        self.name = name
        # Normalized P(E|C)
        self.pec = MCTables.pec
        self.pec_err = MCTables.pec_err
        # Errors on P(E|C), for Vc1
        self.pec_cov_type = "Multinomial"
        # self.pec_cov_type = "Poisson"
        self.PecCov = self.SetPecCovType(self.pec_cov_type)
        # Cause Efficiencies
        self.cEff = MCTables.eff
        self.cEff_inv = safe_inverse(self.cEff)
        # Effective Number of Sim Events
        self.NCmc = MCTables.NCmc
        # Observed Effects Distribution
        self.NEobs = data
        self.NEobs_err = data_err
        # Total Number of Effects
        self.nobs = np.sum(self.NEobs)

        '''Useful Variables for Calculations'''
        # Number of Cause and Effect Bins
        dims = self.pec.shape
        self.cbins = dims[1]
        self.ebins = dims[0]
        # Mixing Matrix
        self.Mij = np.zeros(dims)
        # Flag to Propagate Derivative
        self.ErrorPropFlag = self.SetErrorPropFlag(self.name)
        # Adye Propagating Derivative
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
        #  dcdn = Mij
        dcdn = Mij.copy()
        #  dcdP = ...
        dcdP = np.zeros(self.dcdP.shape)
        f_inv = safe_inverse(f_norm)

        NE_F_R = self.NEobs * f_inv
        # (ti, ec_j + tk) elements
        for ej, ti, tk in product(range(0, ebins), range(0, cbins), range(0, cbins)):
            ec_j = ej*cbins
            dcdP[ti,ec_j+tk] = -NE_F_R[ej] * Mij[ej,ti] * n_c_prev[tk]
        # (ti, ec_j + ti) elements
        for ej, ti in product(range(0, ebins), range(0, cbins)):
            ec_j = ej*cbins
            dcdP[ti,ec_j+ti] += (n_c_prev[ti]*NE_F_R[ej]-n_c[ti])*self.cEff_inv[ti]

        # Adye Propagation Corrections
        if self.ErrorPropFlag and self.counter > 0:
            # Get Previous Derivatives
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
            #  My Version (from Unfolding doc)
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
        # Get Derivative
        dcdn = self.dcdn.copy()
        # Get NObs Covariance
        Vcd = self.getVcd()
        # Set Data Covariance
        Vc0 = dcdn.T.dot(Vcd).dot(dcdn)

        return Vc0

    def getVcPP(self):
        """Get Covariance Matrix of P(E|C), ie from MC
        """
        cbins = self.cbins
        ebins = self.ebins

        CovPP = np.zeros((cbins * ebins, cbins * ebins))

        # Poisson Covariance Matrix
        if self.PecCov == 1:
            for ej, ti in product(range(0, ebins), range(0, cbins)):
                ejc = ej*cbins
                CovPP[ejc+ti,ejc+ti] = self.pec_err[ej,ti]**2
        # Multinomial Covariance Matrix
        elif self.PecCov == 2:

            NC_inv = safe_inverse(self.NCmc)

            for ej, ti in product(range(0, ebins), range(0, cbins)):
                ejc = ej * cbins
                # Don't go looping if zeros are present :)
                if (self.pec[ej,ti] > 0 and NC_inv[ti] > 0):
                    CovPP[ejc+ti,ejc+ti] = NC_inv[ti]*self.pec[ej,ti]*(1-self.pec[ej,ti])
                    cov1 = -NC_inv[ti]*self.pec[ej,ti]
                    for ek in range(ej+1,ebins):
                        ekc = ek*cbins
                        if (ejc+ti == ekc+ti or ek == ej):
                            continue
                        cov = cov1*self.pec[ek,ti]
                        CovPP[ejc+ti,ekc+ti] = cov
                        CovPP[ekc+ti,ejc+ti] = cov

        return CovPP

    def getVc1(self):
        """Get full Vc1 (MC) contribution to cov matrix
        """
        # Get NObs Covariance
        CovPP = self.getVcPP()
        # Get Derivative
        dcdP = self.dcdP
        # Set MC Covariance
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
        # ACM Propagate Errors as Adye
        # DCM Propagate Errors as DAgostini
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
