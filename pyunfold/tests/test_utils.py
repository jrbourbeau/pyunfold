
from __future__ import division, print_function
import numpy as np
import pytest

from pyunfold.Utils import none_to_empty_list, safe_inverse, get_ts


def test_none_to_empty_list_single_input():
    assert none_to_empty_list(None) == []
    assert none_to_empty_list('not empty') == 'not empty'


def test_none_to_empty_list_multi_input():
    a_original, b_original, c_original = 1, None, 'a string'
    a, b, c = none_to_empty_list(a_original, b_original, c_original)

    assert a == a_original
    assert c == c_original
    assert isinstance(b, list) and len(b) == 0


@pytest.mark.parametrize('dtype', [int, float])
def test_safe_inverse(dtype):
    a = np.array([1, 0, 3, 0, 5], dtype=dtype)
    a_inv = safe_inverse(a)
    is_zero = a == 0
    for idx, value in enumerate(a_inv):
        if a[idx] == 0:
            assert value == 0
        else:
            assert value != 0


def test_get_ts_raises():
    with pytest.raises(ValueError) as excinfo:
        get_ts(name='not a valid name')
    assert 'Invalid test statisitc' in str(excinfo.value)
