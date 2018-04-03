
import pandas as pd


class IterativeUnfolder(object):
    """Common base class for iterative unfolder
    """
    def __init__(self, n_c=None, mixer=None, ts_func=None, max_iter=100):

        self.max_iter = max_iter
        # Mixing, Regularizing, Test Statistic Functions
        self.Mix = mixer
        # self.Rglzr = reg_func
        self.ts_func = ts_func
        self.current_n_c = n_c.copy()
        self.counter = 0

    def continue_unfolding(self):
        continue_unfolding = not self.ts_func.pass_tol() and self.counter < self.max_iter
        return continue_unfolding

    def unfold(self):
        """Perform iterative unfolding

        Parameters
        ----------
        self : self

        Returns
        -------
        unfolding_result : pandas.DataFrame
            DataFrame containing the unfolded result for each iteration.
            Each row in unfolding_result corresponds to an iteration.
        """
        unfolding_results = []

        while self.continue_unfolding():
            # Updated unfolded distribution
            # n_c = unfolded_n_c.copy()
            # Mix w/n_c from previous iter
            unfolded_n_c = self.Mix.smear(self.current_n_c)

            # Add mixing result to unfolding_result
            status = {'unfolded': unfolded_n_c,
                      'stat_err': self.Mix.get_stat_err(),
                      'sys_err': self.Mix.get_MC_err()}
            unfolding_results.append(status)

            TS_cur, TS_del, TS_prob = self.ts_func.GetStats(unfolded_n_c,
                                                            self.current_n_c)
            self.current_n_c = unfolded_n_c.copy()
            self.counter += 1

        # Convert unfolding_result dictionary to a pandas DataFrame
        columns = ['sys_err', 'unfolded', 'stat_err']
        unfolding_results = pd.DataFrame.from_records(unfolding_results,
                                                      columns=columns)

        return unfolding_results
