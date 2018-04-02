"""
   Functions to load monte carlo statistic
   data from provided root file.
"""
import numpy as np
import RootReader as rr
from .Utils import safe_inverse, none_to_empty_list


np.set_printoptions(threshold='nan')


class MCTables(object):
    """Base class for loading Monte Carlo tables (Eff, ResponseMatrix)
    """
    def __init__(self, MCFile, BinName=None, EffName="", RespMatrixName="",
                 Stack=False, **kwargs):

        BinName = none_to_empty_list(BinName)

        # Strings Naming Things
        self.MCFile = MCFile
        self.BinName = BinName
        # Flag to Stack
        self.StackFlag = Stack
        self.nStack = len(BinName)
        # Empty Objects for Stats and Axes
        self.NCmc = None
        self.eff = None
        self.pec = None
        self.pec_err = None
        self.Caxis = [[] for i in range(self.nStack)]
        self.Cedges = [[]for i in range(self.nStack)]
        self.Eaxis = []
        self.Eedges = []
        self.Clabel = ''
        self.Elabel = ''

        # Load Objects if Names Provided
        self.NCLoaded = False
        self.EffLoaded = False
        if EffName:
            self.load_eff_table(EffName)
        self.RespMLoaded = False
        if RespMatrixName:
            self.load_response_matrix(RespMatrixName)

    def print_load_message(self, name, hist_name):
        load_msg = ("\n***{} table, {}, already loaded! Ignoring "
                    "load request.***\n".format(name, hist_name))
        print(load_msg)

    def check_array(self, inarray, att_axis, name):
        # First check if attribute axis is empty
        err_mess = "\n*** {} array is empty! Exiting... ***\n".format(name)
        assert len(inarray) != 0, err_mess
        if att_axis == []:
            out = inarray.copy()
        # If not, ensure that inarray matches atrribute x-axis already loaded.
        else:
            err_mess = ('\n*** Mismatch between two X-axes in MC file {}. '
                        'Exiting... ***\n'.format(self.MCFile))
            assert np.array_equal(inarray, att_axis), err_mess
            out = att_axis.copy()
        return out

    def GetNCmc(self):
        assert self.NCLoaded, "\n*** MC NC table not loaded. Exiting... ***\n"
        return self.NCmc

    def GetEff(self):
        assert self.EffLoaded, "\n*** Effective Area table not loaded. Exiting... ***\n"
        return self.eff

    def GetPEC(self):
        assert self.RespMLoaded, "\n*** Response Matrix not loaded. Exiting... ***\n"
        return self.pec

    def GetPECError(self):
        assert self.RespMLoaded, "\n*** Response Matrix not loaded. Exiting... ***\n"
        return self.pec_err

    def GetCauseLabel(self):
        return self.Clabel

    def GetEffectLabel(self):
        return self.Elabel

    # Get the axis centers and edges for causes
    def GetCauseAxis(self, index):
        err_msg = "\n*** Cause axis not loaded. Exiting... ***\n"
        assert (not (self.Caxis[index] == [] or self.Cedges[index] == [])), err_msg
        return self.Caxis[index], self.Cedges[index]

    # Get the axis centers and edges for effects
    def GetEffectAxis(self):
        err_msg = "\n*** Effect axis not loaded. Exiting... ***\n"
        assert (not (self.Eaxis == [] or self.Eedges == [])), err_msg
        return self.Eaxis, self.Eedges

    # Get the Effective Area
    def load_eff_table(self, name):
        tabname = "Eff"
        if not self.EffLoaded:
            eff_tabs = []
            NCmc_tabs = []
            for ibin in range(self.nStack):
                axis, edges, eff, eff_err = rr.get1d(self.MCFile, name, self.BinName[ibin])
                Clabel, Efflabel, Title = rr.get_labels(self.MCFile, name, self.BinName[ibin])

                eff_tabs.append(eff)
                eff_err_inv = safe_inverse(eff_err)
                NCmc_tabs.append((eff * eff_err_inv)**2)

                self.Caxis[ibin] = self.check_array(axis,
                                                    self.Caxis[ibin],
                                                    tabname + " cause axis")
                self.Cedges[ibin] = self.check_array(edges,
                                                     self.Cedges[ibin],
                                                     tabname + " cause xedges")
                self.EffLoaded = True
                self.NCLoaded = True
            if self.nStack == 1:
                self.eff = eff_tabs[0]
                self.NCmc = NCmc_tabs[0]
            else:
                self.eff = np.hstack(eff_tabs)
                self.NCmc = np.hstack(NCmc_tabs)
        else:
            self.print_load_message(tabname, name)

    # Get the Migration Matrix table
    def load_response_matrix(self, name):
        tabname = "Resp Matrix"
        if not self.RespMLoaded:
            pec_tabs = []
            pec_err_tabs = []
            for ibin in range(self.nStack):
                axes, edges, pec, pec_err = rr.get2d(self.MCFile, name, self.BinName[ibin])
                Clabel, Elabel, Title = rr.get_labels(self.MCFile, name, self.BinName[ibin])

                self.Clabel = Clabel
                self.Elabel = Elabel

                Caxis = axes[0].copy()
                Eaxis = axes[1].copy()
                Cedges = edges[0].copy()
                Eedges = edges[1].copy()

                pec_tabs.append(pec)
                pec_err_tabs.append(pec_err)

                self.Caxis[ibin] = self.check_array(Caxis,
                                                    self.Caxis[ibin],
                                                    tabname + " cause axis")
                self.Cedges[ibin] = self.check_array(Cedges,
                                                     self.Cedges[ibin],
                                                     tabname + " cause edges")
                self.Eaxis = self.check_array(Eaxis, self.Eaxis, tabname + " effects axis")
                self.Eedges = self.check_array(Eedges, self.Eedges, tabname + " effects edges")
                self.RespMLoaded = True
            if (self.nStack == 1):
                self.pec = pec_tabs[0]
                self.pec_err = pec_err_tabs[0]
            else:
                self.pec = np.hstack(pec_tabs)
                self.pec_err = np.hstack(pec_err_tabs)
        else:
            self.print_load_message(tabname, name)
