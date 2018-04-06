
from __future__ import division, print_function
import pytest

from pyunfold.teststat import get_ts


def test_get_ts_raises():
    with pytest.raises(ValueError) as excinfo:
        get_ts(name='not a valid name')
    assert 'Invalid test statisitc' in str(excinfo.value)
