
from __future__ import division, print_function
from functools import wraps
import inspect
import numpy as np


def assert_kwargs_not_none(*inputs):
    """Decorator to wrap functions that require non-None kwargs

    Parameters
    ----------
    inputs : str
        Name of keyword arguments to check for non-None values.

    Raises
    ------
    ValueError
        If one of the ``inputs`` keyword arguments is None.
    """
    def real_decorator(func):
        # Want to build list of keyword arguments from the func signature
        argspec = inspect.getargspec(func)
        sig_kwarg_keys = argspec.args[-len(argspec.defaults):]
        sig_kwarg_values = argspec.defaults
        sig_kwargs = dict(zip(sig_kwarg_keys, sig_kwarg_values))

        @wraps(func)
        def wrapper(*args, **kwargs):
            for input_ in inputs:
                err_msg = 'Input keyword argument "{}" must not be None'.format(input_)
                try:
                    if kwargs[input_] is None:
                        raise ValueError(err_msg)
                except KeyError:
                    if sig_kwargs[input_] is None:
                        raise ValueError(err_msg)
            return func(*args, **kwargs)
        return wrapper
    return real_decorator


def assert_same_shape(*arrays):
    """Checks that each input array_like objects are the same shape
    """
    arrays = cast_to_array(*arrays)
    shapes = [array.shape for array in arrays]
    unique_shapes = set(shapes)
    if not len(unique_shapes) == 1:
        raise ValueError('Multiple shapes found: {}'.format(unique_shapes))


def cast_to_array(*arrays):
    """Casts input arrays to numpy.ndarray objects

    Note that no copy is made if an input array is already a numpy.ndarray.

    Parameters
    ----------
    arrays : array_like
        Input array_like objects to be cast to numpy arrays.

    Returns
    -------
    output : list
        List of casted numpy arrays.

    Examples
    --------
    >>> import numpy as np
    >>> a_original = [1, 2, 3]
    >>> b_original = np.array([4.5, 2.1, 900])
    >>> a, b = cast_to_array(a_original, b_original)
    >>> a
    array([1, 2, 3])
    >>> b
    array([  4.5,   2.1, 900. ])
    >>> b is b_original
    True
    """
    if len(arrays) == 1:
        output = np.asarray(arrays[0])
    else:
        output = map(np.asarray, arrays)
    return output


def none_to_empty_list(*args):
    """Replaces None inputs with an empty list

    Examples
    --------
    Single input case

    >>> none_to_empty_list(None)
    []

    Multiple input case

    >>> a, b, c = None, 'woo', 34
    >>> none_to_empty_list(a, b, c)
    [[], 'woo', 34]
    """
    outputs = []
    for arg in args:
        outputs.append(arg if arg is not None else [])
    if len(outputs) == 1:
        return outputs[0]
    else:
        return outputs


def safe_inverse(x):
    """Safely inverts the elements in x

    Parameters
    ----------
    x : array_like
        Input array to take the inverse of (i.e. 1 / x).

    Returns
    -------
    inv : numpy.ndarray
        Inverse of input array with inf set to zero.

    Examples
    --------
    >>> a = [1, 2, 3, 0, 4]
    >>> safe_inverse(a)
    array([1.        , 0.5       , 0.33333333, 0.        , 0.25      ])
    """
    x = np.asarray(x)
    is_zero = x == 0
    with np.errstate(divide='ignore'):
        inv = 1 / x
    inv[is_zero] = 0

    return inv
