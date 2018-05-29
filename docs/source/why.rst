.. _features:

:github_url: https://github.com/jrbourbeau/pyunfold

*************
Why PyUnfold?
*************

PyUnfold was built out of a need for an easy to use, yet extensible statistical
deconvolution package. It is built on top of the Python scientific
computing stack, bringing iterative unfolding methods fully into the
Python ecosystem.

-----
Goals
-----

PyUnfold has been designed to be easy to use for first-time users, flexible
enough to perform real-world analyses, and allows for seamlessly testing the
robustness of results. Users can perform an unfolding with just a few lines of
code

.. code-block:: python

    from pyunfold import iterative_unfold

    unfolded_result = iterative_unfold(data=data,
                                       data_err=data_err,
                                       response=response,
                                       response_err=response_err,
                                       efficiencies=efficiencies,
                                       efficiencies_err=efficiencies_err)


PyUnfold provides users:

- A tool for the analysis of measurement effects on physical distributions
- A straightforward API that's suitable for many experimental applications
- The ability to easily incorporate both statistical and systematic uncertainties


-------------
Main Features
-------------

The features unique to PyUnfold include:

- Built on top of the Python scientific computing stack (i.e. NumPy, SciPy, pandas)
- Support for custom, user defined initial prior distributions
- Optional regularization via tunable univariate splines
- Adjustable test statistic convergence criterion comparing unfolded distributions between iterations
- Support for multi-dimensional unfolding


---------
Successes
---------

PyUnfold has been successfully used in several contexts, including:

- Cosmic-ray energy spectrum measurement made by the `HAWC observatory <https://www.hawc-observatory.org/>`_ [1]_.
- Cosmic-ray composition analysis using the `IceCube South Pole Neutrino Observatory <https://icecube.wisc.edu/>`_.


.. [1] Alfaro, R. and others. 2017. "All-particle cosmic ray energy spectrum measured by the HAWC experiment from 10 to 500 TeV." *Phys. Rev. D* 96 (12):122001. `<https://doi.org/10.1103/PhysRevD.96.122001>`_.
