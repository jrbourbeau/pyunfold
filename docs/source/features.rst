.. _features:

:github_url: https://github.com/jrbourbeau/pyunfold

*************
Why PyUnfold?
*************

PyUnfold was built out of a need to facilitate statistical deconvolution in a concise
and convenient, yet extensible package.
Born within the realm of experimental high-energy physics (HEP), PyUnfold removes
dependencies on HEP-specific software, instead having its structure rest on the 
the Python scientific computing stack.
This not only simplifies unfolding for the wide range of physics analyses, but also
permits its implementation beyond the scope of the HEP community.


-----
Goals
-----

PyUnfold has been designed to be both easy to use for first-time users as well as flexible 
enough for fine-tuning an analysis and seamlessly testing the robustness of results.
PyUnfold provides users:

- A tool for the analysis of measurement effects on physical distributions
- A straightforward API that's suitable for many experimental applications
- The ability to easily incorporate both statistical and systematic uncertainties


-------------
Main Features
-------------

The features unique to PyUnfold are the following:

- Dependency on Python scientific computing stack (NumPy, SciPy, Pandas)

- Flexibility to test results based on user provided input:

  - Custom, user defined initial prior distributions

  - Optional regularization via tunable univariate splines

  - Adjustable test statistic convergence criterion comparing distributions between iterations

- Multidimensional unfolding in cause subsets, which are regularized in their respective groups


PyUnfold has been used successfully in astroparticle physics measurements by both
the `HAWC <https://www.hawc-observatory.org/>`_ and `IceCube <https://icecube.wisc.edu/>`_
observatories.
