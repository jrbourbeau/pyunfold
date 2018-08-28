
import pytest

from pyunfold.teststat import get_ts, TEST_STATISTICS


@pytest.mark.parametrize('name', TEST_STATISTICS.keys())
def test_get_ts(name):
    assert get_ts(name) in TEST_STATISTICS.values()


def test_get_ts_raises():
    name = 'not a valid ts'
    with pytest.raises(ValueError) as excinfo:
        get_ts(name)

    expected_msg = ('Invalid test statistic, {}, entered. Must be '
                    'in {}'.format(name, TEST_STATISTICS.keys()))
    assert expected_msg == str(excinfo.value)


@pytest.mark.parametrize('ts', TEST_STATISTICS.keys())
def test_ts_calc(ts, example_dataset):
    # Regression test for issue #92
    ts_obj = get_ts(ts)
    ts_func = ts_obj(tol=0.01,
                     num_causes=len(example_dataset.data),
                     TestRange=[0, 1e2],
                     verbose=False)
    ts_func.calc(example_dataset.data, example_dataset.data + 1)
