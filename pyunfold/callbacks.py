
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

    Smooths the unfolded distribution at each iteration using
    ``UnivariateSpline`` from ``scipy.interpolate``. For more information about
    ``UnivariateSpline``, see the
    `UnivariateSpline API documentation
    <https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html>`_.

    Parameters
    ----------
    degree : int, optional
        Degree of the smoothing spline. Must be <= 5 (default is 3, a cubic
        spline).
    smooth : float or None, optional
        Positive smoothing factor used to choose the number of knots. If 0,
        spline will interpolate through all data points (default is None).
    groups : array_like, optional
        Group labels for each cause bin. If groups are specified, then each
        cause group will be regularized independently (default is None).

    Notes
    -----
    The number of causes must be larger than the spline ``degree``.

    Examples
    --------
    Specify the spline degree and smoothing factor:

    >>> from pyunfold.callbacks import SplineRegularizer
    >>> reg = SplineRegularizer(degree=3, smooth=1.25)

    Different cause groups are also supported. For instance, in a problem with
    seven cause bins, if the first three cause bins belong to their own group,
    the next two cause bins belong to another group, and the last two cause
    bins belong to yet another group, an array can be constructed that
    identifies the group each cause bin belongs to. E.g.

    >>> groups = [0, 0, 0, 1, 1, 2, 2]
    >>> reg = SplineRegularizer(degree=3, smooth=1.25, groups=groups)

    If provided with a ``groups`` parameter, ``SplineRegularizer`` will
    regularize the unfolded distribution for each group independently.
    """
    def __init__(self, degree=3, smooth=None, groups=None):
        super(SplineRegularizer, self).__init__()
        self.degree = degree
        self.smooth = smooth
        self.groups = np.asarray(groups) if groups is not None else None

    def on_iteration_end(self, iteration, params):
        y = params['unfolded']
        x = np.arange(len(y), dtype=float)
        if self.groups is None:
            spline = UnivariateSpline(x, y, k=self.degree, s=self.smooth)
            fitted_unfolded = spline(x)
        else:
            # Check that a group is specified for each cause bin
            if len(self.groups) != len(y):
                err_msg = ('Invalid groups array. There should be an entry '
                           'for each cause bin. However, got len(groups)={} '
                           'while there are {} cause bins.'.format(len(self.groups), len(y)))
                raise ValueError(err_msg)
            fitted_unfolded = np.empty(len(y))
            group_ids = np.unique(self.groups)
            for group in group_ids:
                group_mask = self.groups == group
                x_group = x[group_mask]
                y_group = y[group_mask]
                spline_group = UnivariateSpline(x_group, y_group, k=self.degree, s=self.smooth)
                fitted_unfolded_group = spline_group(x_group)
                fitted_unfolded[group_mask] = fitted_unfolded_group

        return fitted_unfolded


def validate_callbacks(callbacks):
    if callbacks is None:
        callbacks = []
    elif isinstance(callbacks, Callback):
        callbacks = [callbacks]
    else:
        if not all([isinstance(c, Callback) for c in callbacks]):
            invalid_callbacks = [c for c in callbacks if not isinstance(c, Callback)]
            raise TypeError('Found non-callback object in callbacks: {}'.format(invalid_callbacks))

    return callbacks


def extract_regularizer(callbacks):
    callbacks = validate_callbacks(callbacks)
    regularizers = [c for c in callbacks if isinstance(c, Regularizer)]
    if len(regularizers) > 1:
        raise NotImplementedError('Multiple regularizer callbacks where provided.')
    regularizer = regularizers[0] if len(regularizers) == 1 else None

    return regularizer
