#!/usr/bin/env python
"""
   Covariance Matrix Calculator for D'Agostini
   unfolding method, with modifications to
   error propagation based on method request
   from user input.

.. codeauthor: Zig Hampel
"""

__version__ = "$Id"


try:
    import numpy as np
    from Utils import safe_inverse
except ImportError as e:
    print e
    raise ImportError

class CovarianceMatrix:
    'Covariance Base Class'
    def __init__(self, name, MCTables, EffectsDist):

        '''From Initialization Inputs'''
        self.name = name
        # Normalized P(E|C)
        self.pec = MCTables.GetPEC()
        self.pec_err = MCTables.GetPECError()
        # Errors on P(E|C), for Vc1
        self.pec_cov_type = "Multinomial"
        #self.pec_cov_type = "Poisson"
        self.PecCov = self.SetPecCovType(self.pec_cov_type)
        # Cause Efficiencies
        self.cEff = MCTables.GetEff()
        self.cEff_inv = safe_inverse(self.cEff)
        # Effective Number of Sim Events
        self.NCmc = MCTables.GetNCmc()
        # Observed Effects Distribution
        self.NEobs = EffectsDist.getData()
        self.NEobs_err = EffectsDist.getError()
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
        self.dcdP = np.zeros((self.cbins,self.cbins*self.ebins))
        # Counter for number of iterations
        self.counter = 0

    """
        All Covariance Matrix Code either transcribed from Adye's RooUnfoldBayes.cxx
        (http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html#dev),
        or based on derivations from Unfolding HAWC reference 'Unfolding Uncertainties'
        (http://private.hawc-observatory.org/hawc.umd.edu/internal/doc.php?id=2329).
    """
    
    '''Set the Current State of dcdn and dcdP'''
    def setCurrent_State(self, Mij, f_norm, n_c, n_c_prev):
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

        NE_F_R = self.NEobs*f_inv
        for ej in xrange(0,ebins):
            ne_f_r = NE_F_R[ej]
            # Combined counting index
            ec_j = ej*cbins
            for ti in xrange(0,cbins):
                b = -ne_f_r * Mij[ej,ti]
                for tk in xrange(0,cbins):
                    dcdP[ti,ec_j+tk] = b*n_c_prev[tk]
                val = n_c_prev[ti]*ne_f_r-n_c[ti]
                dcdP[ti,ec_j+ti] += (n_c_prev[ti]*ne_f_r-n_c[ti])*self.cEff_inv[ti]

        # Adye Propagation Corrections
        if (self.ErrorPropFlag and self.counter > 0):
            # Get Previous Derivatives
            dcdn_prev = self.dcdn.copy()
            dcdP_prev = self.dcdP.copy()

            n_c_prev_inv = safe_inverse(n_c_prev)
            # Ratio of updated n_c to n_c_prev
            nc_r = n_c*n_c_prev_inv
            # Efficiency ratio of n_c_prev
            e_r = self.cEff*n_c_prev_inv

            # Calculate extra dcdn terms
            M1 = dcdn_prev.copy()
            M2 = Mij.copy()
            for tj in xrange(0,cbins):
                M1[:,tj] *= nc_r[tj]
                M2[:,tj] *= -e_r[tj]
            for ej in xrange(0,ebins):
                M2[ej,:] *= self.NEobs[ej]
            M3 = np.dot(M2.T,dcdn_prev)
            dcdn += np.dot(Mij,M3)
            dcdn += M1

            # Calculate extra dcdP terms
            #  My Version (from Unfolding HAWC-doc)
            A = Mij.copy()
            B = Mij.copy()
            for ej in xrange(0,ebins):
                A[ej,:] *= self.NEobs[ej]
            for ti in xrange(0,cbins):
                B[:,ti] *= e_r[ti]
            C = np.dot(A.T,B)
            dcdP_Upd = np.dot(C,dcdP_prev)
            nec = ebins*cbins
            for tj in xrange(0,cbins):
                r = nc_r[tj]
                for jk in xrange(0,nec):
                    dcdP[tj,jk] += r*dcdP_prev[tj,jk] - dcdP_Upd[tj,jk]

        # Set current derivative matrices
        self.dcdn = dcdn.copy()
        self.dcdP = dcdP.copy()
        # On to next iteration
        self.counter += 1
    
    ''' Get Covariance Matrix of N(E), ie from Observed Effects '''
    def getVcd(self):

        ebins = self.ebins
        Vcd = np.zeros((ebins,ebins))

        for ei in xrange(0,ebins):
            Vcd[ei,ei] = self.NEobs_err[ei]**2

        return Vcd
    
    ''' Get full Vc0 (data) contribution to cov matrix '''
    def getVc0(self):

        ebins = self.ebins
        cbins = self.cbins

        # Get Derivative
        dcdn = self.dcdn.copy()
        # Get NObs Covariance
        Vcd = self.getVcd()
        # Set Data Covariance 
        Vc0 = dcdn.T.dot(Vcd).dot(dcdn)
        return Vc0
                
    ''' Get Covariance Matrix of P(E|C), ie from MC'''
    def getVcPP(self):
        cbins = self.cbins
        ebins = self.ebins

        CovPP = np.zeros((cbins*ebins,cbins*ebins))

        # Poisson Covariance Matrix
        if (self.PecCov == 1):
            for ej in xrange(0,ebins):
                ejc = ej*cbins
                for ti in xrange(0,cbins):
                    CovPP[ejc+ti,ejc+ti] = self.pec_err[ej,ti]**2
        # Multinomial Covariance Matrix
        elif (self.PecCov == 2):

            NC_inv = safe_inverse(self.NCmc)

            for ej in xrange(0,ebins):
                ejc = ej*cbins
                for ti in xrange(0,cbins):
                    # Don't go looping if zeros are present :)
                    if (self.pec[ej,ti] > 0 and NC_inv[ti] > 0):
                        CovPP[ejc+ti,ejc+ti] = NC_inv[ti]*self.pec[ej,ti]*(1-self.pec[ej,ti])
                        cov1 = -NC_inv[ti]*self.pec[ej,ti]
                        for ek in xrange(ej+1,ebins):
                            ekc = ek*cbins
                            if (ejc+ti == ekc+ti or ek == ej):
                                continue
                            cov = cov1*self.pec[ek,ti]
                            CovPP[ejc+ti,ekc+ti] = cov
                            CovPP[ekc+ti,ejc+ti] = cov

        return CovPP

    ''' Get full Vc1 (MC) contribution to cov matrix '''
    def getVc1(self):
        # Get NObs Covariance
        CovPP = self.getVcPP()
        # Get Derivative
        dcdP = self.dcdP
        # Set MC Covariance
        Vc1 = dcdP.dot(CovPP).dot(dcdP.T)
        return Vc1
    
    ''' Get full covariance matrix '''
    def getCov(self):
        # Get Contributions from Vc0 and Vc1
        Vc0 = self.getVc0()
        Vc1 = self.getVc1()
        # Full Covariance Matrix
        Vc = Vc0+Vc1
        return Vc

    # Keys and Chooser Funtions for Covariance Matrix Options

    ''' Key to choose whether to propagate errors at each iteration '''
    ErrorPropOptionsKey = {
        "ACM" : True, #Propagate Errors as Adye
        "DCM" : False #Propagate Errors as DAgostini
        }
    
    def SetErrorPropFlag(self,func=None):
        if func in self.ErrorPropOptionsKey:
            return self.ErrorPropOptionsKey[func]
        else:
            print("==================================")
            print("Covariance calculation method not found.")
            print("Cov matrix options are:")
            for keys,values in self.ErrorPropOptionsKey.items():
                print("  "+keys+":\t"+values.__doc__)
            print("Exiting...")
            import sys
            sys.exit(0)
    
    ''' Key for Selecting PEC Covariance Matrix Type ''' 
    PecCovOptionsKey = {
        "Poisson" : 1,
        "Multinomial" : 2
        }
    
    def SetPecCovType(self,func=None):
        if func in self.PecCovOptionsKey:
            return self.PecCovOptionsKey[func]
        else:
            print("==================================")
            print("Covariance calculation method not found.")
            print("Cov matrix options are:")
            for keys,values in self.PecCovOptionsKey.items():
                print("  "+keys+":\t"+values.__doc__)
            print("Exiting...")
            import sys
            sys.exit(0)
