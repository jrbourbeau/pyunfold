
import pytest

from pyunfold.teststat import get_ts, TEST_STATISTICS


@pytest.mark.parametrize('name', ['ks',
                                  'chi2',
                                  'pf',
                                  'rmd'])
def test_get_ts(name):
    assert name in TEST_STATISTICS.keys()
    assert get_ts(name) in TEST_STATISTICS.values()


def test_get_ts_raises():
    name = 'not a valid ts'
    with pytest.raises(ValueError) as excinfo:
        get_ts(name)

    expected_msg = ('Invalid test statisitc, {}, entered. Must be '
                    'in {}'.format(name, TEST_STATISTICS.keys()))
    assert expected_msg == str(excinfo.value)
