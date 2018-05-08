
from __future__ import division, print_function
import numpy as np
import pytest

from pyunfold.utils import (none_to_empty_list, safe_inverse, cast_to_array,
                            assert_same_shape)


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
    for idx, value in enumerate(a_inv):
        if a[idx] == 0:
            assert value == 0
        else:
            assert value != 0


def test_cast_to_array_multi_input():
    a_original = [1, 2, 3]
    b_original = np.array([4, 5, 6])

    a, b = cast_to_array(a_original, b_original)

    for output, input_ in zip([a, b], [a_original, b_original]):
        assert isinstance(output, np.ndarray)
        np.testing.assert_allclose(output, input_)


def test_cast_to_array_single_input():
    a_original = [1, 2, 3]
    a = cast_to_array(a_original)
    assert isinstance(a, np.ndarray)
    np.testing.assert_allclose(a_original, a)


def test_cast_to_array_no_copy():
    # Check that input numpy arrays are not copied
    a_original = np.array([1, 2, 3])
    a = cast_to_array(a_original)
    assert a is a_original


def test_assert_same_shape():
    a = [1, 2, 3]
    b = [4, 5, 6]
    assert_same_shape(a, b)


def test_assert_same_shape_raises_1d():
    a = [1, 2, 3]
    b = [4, 5, 6, 7]
    with pytest.raises(ValueError):
        assert_same_shape(a, b)


def test_assert_same_shape_raises_2d():
    a = [[1, 2, 3], [4, 5, 6]]
    b = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    with pytest.raises(ValueError):
        assert_same_shape(a, b)
