
import pyunfold


def test_version_exists():
    assert hasattr(pyunfold, '__version__')


def test_iterative_unfold_exists():
    assert hasattr(pyunfold, 'iterative_unfold')
