
import sys
import numpy as np
from scipy.interpolate import UnivariateSpline


class Callback(object):
    """Callback base class
    """
    def __init__(self):
        pass

    def on_unfolding_begin(self):
        pass

    def on_unfolding_end(self):
        pass

    def on_iteration_begin(self, iteration, params):
        pass

    def on_iteration_end(self, iteration, params):
        pass


class Regularizer(object):
    """Regularizer callback base class
    """
    def __init__(self):
        pass


class Logger(Callback):
    """Logger callback

    Writes test statistic information for each iteration to sys.stdout.
    """
    def __init__(self):
        super(Callback, self).__init__()
        pass

    def on_iteration_end(self, iteration, params):
        """Writes to sys.stdout

        Parameters
        ----------
        iteration : int
            Unfolding iteration (i.e. iteration=1 after first unfolding
            iteration, etc.).
        params : dict
            Dictionary containing key value pairs for the current test
            statistic value (``'ts_iter'``) and the final test statistic stopping
            condition (``'ts_stopping'``).
        """
        output = ('Iteration {}: ts = {:0.4f}, ts_stopping ='
                  ' {}\n'.format(iteration,
                                 params['ts_iter'],
                                 params['ts_stopping']))
        sys.stdout.write(output)


class SplineRegularizer(Callback, Regularizer):
    """Spline regularization callback

    Fits scipy.interpolate.UnivariateSpline to unfolded distribution at each
    iteration.
    """
    def __init__(self, degree=3, smooth=None):
        super(SplineRegularizer, self).__init__()
        self.degree = degree
        self.smooth = smooth

    def on_iteration_end(self, iteration, params):
        y = params['unfolded']
        x = np.arange(len(y), dtype=float)
        spline = UnivariateSpline(x, y, k=self.degree, s=self.smooth)
        fitted_unfolded = spline(x)

        return fitted_unfolded
