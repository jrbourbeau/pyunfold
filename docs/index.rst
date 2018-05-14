.. PyUnfold documentation master file, created by
   sphinx-quickstart on Thu Mar 29 14:30:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:github_url: https://github.com/jrbourbeau/pyunfold

PyUnfold
========

.. image:: https://travis-ci.org/jrbourbeau/pyunfold.svg?branch=master
    :target: https://travis-ci.org/jrbourbeau/pyunfold

.. image:: https://ci.appveyor.com/api/projects/status/wphmmposuctye5ye/branch/master?svg=true
    :target: https://ci.appveyor.com/project/jrbourbeau/pyunfold/branch/master

.. image:: https://codecov.io/gh/jrbourbeau/pyunfold/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/jrbourbeau/pyunfold

.. image:: https://img.shields.io/pypi/v/pyunfold.svg
    :target: https://pypi.org/project/PyUnfold/

.. image:: https://img.shields.io/pypi/pyversions/pyunfold.svg
    :target: https://pypi.org/project/PyUnfold/

.. image:: https://img.shields.io/pypi/l/pyunfold.svg
    :target: https://pypi.org/project/PyUnfold/

PyUnfold is a Python package for implementing the D'Agostini iterative unfolding algorithm [1]_.

Quickstart
----------

Implementing iterative unfolding with PyUnfold only takes a few lines of code.

.. code-block:: python

    from pyunfold import iterative_unfold

    # Counts distributions
    data = [100, 150]
    data_err = [10, 12.2]

    # Response matrix
    response = [[0.9, 0.1],
                [0.1, 0.9]]
    response_err = [[0.01, 0.01],
                    [0.01, 0.01]]

    # Detection efficiencies
    efficiencies = [1, 1]
    efficiencies_err = [0.01, 0.01]

    # Perform iterative unfolding
    unfolded_result = iterative_unfold(data, data_err,
                                       response, response_err,
                                       efficiencies, efficiencies_err)

The returned unfolded result is a dictionary containing:

- ``unfolded``: Unfolded cause distribution
- ``stat_err``: Statistical (Poisson) errors
- ``sys_err``: Systematic errors associated with limited statistics in the response matrix

.. code-block:: python

    {'unfolded': array([ 94.48002622, 155.51997378]),
     'sys_err': array([0.66204237, 0.6620424]),
     'stat_err': array([11.2351567, 13.75617997])}



.. toctree::
    :maxdepth: 1
    :caption: Getting Started

    overview
    installation
    Tutorial <tutorial.ipynb>

.. toctree::
    :maxdepth: 1
    :caption: User Guide

    callbacks
    api
    advanced
    mathematics
    changelog
    contributing

.. toctree::
   :maxdepth: 1
   :caption: Useful links
   :hidden:

   PyUnfold @ GitHub <https://github.com/jrbourbeau/pyunfold>
   PyUnfold @ PyPI <https://pypi.org/project/pyunfold/>
   Issue Tracker <https://github.com/jrbourbeau/pyunfold/issues>


.. [1] G. D'Agostini, “A Multidimensional unfolding method based on Bayes' theorem”, Nucl. Instrum. Meth. A **362** (1995) 487.
