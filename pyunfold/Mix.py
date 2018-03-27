#!/usr/bin/env python
"""
   D'Agostini mixing routine.
   After MC and observed distribution
   are provided, can perform a single
   mixing, given a prior pdf.

.. codeauthor: Zig Hampel
"""

__version__ = "$Id"

try:
    import sys
    import numpy as np
    from CovMatrix import CovarianceMatrix as CovMatrix
    from Utils import safe_inverse
except ImportError as e:
    print e
    raise ImportError


class Mixer:
    'DAgostini Bayesian Mixer Class'
    def __init__(self, name, ErrorType="", MCTables=[], EffectsDist=[], **kwargs):
        
        '''From Initialization Inputs'''
        self.name = name
        # Normalized P(E|C)
        self.pec = MCTables.GetPEC()
        # Cause Efficiencies
        self.cEff = MCTables.GetEff()
        # Observed Effects Distribution
        self.NEobs = EffectsDist.getData()
        # Total Number of Effects
        self.nobs = np.sum(self.NEobs)
        # Error Calculation Type
        self.ErrorType = ErrorType

        '''Useful Variables for Calculations'''
        # Number of Cause and Effect Bins
        dims = self.pec.shape
        self.cbins = dims[1]
        self.ebins = dims[0]
        # Flag for Correctly Sized Arrays
        self.DimFlag = False
        # Inverse of Cause Efficiencies
        self.cEff_inv = safe_inverse(self.cEff)
        # Mixing Matrix
        self.Mij = np.zeros(dims)

        '''Covariance Matrix'''
        self.Cov = CovMatrix(ErrorType, MCTables, EffectsDist)
    
    
    ''' Test for Equal Length Arrays '''
    def CheckDims(self, p_c):
        if ( len(self.NEobs) != self.ebins ):
            print("Size of observed effects array, len(NEobs), not equal \
                   to mixing matrix effects dimension, pec.shape[1].")
            print("%i != %i" %(len(self.NEobs),self.ebins))
            import sys
            sys.exit(0)
        elif ( len(self.cEff) != self.cbins ):
            print self.cEff
            print("Size of effective area array, len(cEff), not equal \
                   to mixing matrix causes dimension, pec.shape[0].")
            print("%i != %i" %(len(self.cEff),self.cbins))
            import sys
            sys.exit(0)
        # Ensure Prior is Properly Sized
        elif ( len(p_c) != self.cbins ):
            print("Size of prior array, len(p_c), not equal \
                   to mixing matrix causes dimension, pec.shape[1].")
            print("%i != %i" %(len(p_c),self.cbins))
            import sys
            sys.exit(0)
        else:
            self.DimFlag = True
        

    '''Useful Variables for Calculations'''

    # Get the Covariance Matrix
    def getCov(self):
        cvm = self.Cov.getCov()
        return cvm
    
    # Get Statistical Errors
    def getStatErr(self):
        cvm = self.Cov.getVc0()
        err = np.sqrt(cvm.diagonal())
        return err

    # Get MC (Systematic) Errors 
    def getMCErr(self):
        cvm = self.Cov.getVc1()
        err = np.sqrt(cvm.diagonal())
        return err
    
    '''Smear Calculation - the Bayesian Mixer!'''
    # Only needs input prior, n_c, to unfold
    def Smear(self, n_c):
        # Test Input Array Sizes
        if (self.DimFlag == False):
            self.CheckDims(n_c)
        
        ebins = self.ebins
        cbins = self.cbins

        # Bayesian Normalization Term (denominator)
        f_norm = np.dot(self.pec,n_c)
        f_inv = safe_inverse(f_norm)
    
        # Unfolding (Mij) Matrix at current step
        Mij = np.zeros(self.Mij.shape)

        n_c_eff = n_c*self.cEff_inv
        for ti in xrange(0,cbins):
            for ej in xrange(0,ebins):
                Mij[ej,ti] = self.pec[ej,ti]*n_c_eff[ti]*f_inv[ej]
    
        # Estimate cause distribution via Mij
        n_c_update = np.dot(self.NEobs,Mij)

        # The status quo
        self.Mij = Mij.copy()
        self.Cov.setCurrent_State(self.Mij,f_norm,n_c_update,n_c)
       
        return n_c_update
