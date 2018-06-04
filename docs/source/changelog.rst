.. _changelog:

:github_url: https://github.com/jrbourbeau/pyunfold

*************
Release Notes
*************

Version 0.5.dev0 (TBD)
----------------------

**New Features**:

-


**Changes**:

-


**Bug Fixes**:

-


Version 0.4 (2018-06-03)
------------------------

**Changes**:

- Adds check to ``iterative_unfold`` for negative values in input
  distributions. An error is raised if negative efficiencies, responses, or
  observed data distributions are input.
  (See `PR #73 <https://github.com/jrbourbeau/pyunfold/pull/73>`_)
- Adds several documentation updates in response to JOSS reviewer comments
  (See
  `PR #83 <https://github.com/jrbourbeau/pyunfold/pull/83>`_,
  `PR #82 <https://github.com/jrbourbeau/pyunfold/pull/82>`_,
  `PR #77 <https://github.com/jrbourbeau/pyunfold/pull/77>`_,
  `PR #76 <https://github.com/jrbourbeau/pyunfold/pull/76>`_,
  `PR #71 <https://github.com/jrbourbeau/pyunfold/pull/71>`_,
  `PR #70 <https://github.com/jrbourbeau/pyunfold/pull/70>`_, and
  `PR #69 <https://github.com/jrbourbeau/pyunfold/pull/69>`_)


**Bug Fixes**:

- Fixes environment used for running the "Tutorial" and "Advanced Techniques"
  example notebooks on binder. (See `PR #68 <https://github.com/jrbourbeau/pyunfold/pull/68>`_)


Version 0.3 (2018-05-15)
------------------------

**New Features**:

- Adds cause bin group support to ``SplineRegularizer`` for independent
  regularization. (See `PR #47 <https://github.com/jrbourbeau/pyunfold/pull/47>`_)
- Adds ``cov_type`` parameter to the ``iterative_unfold`` function. This adds
  the option to choose between using the Multinomial or Poisson form of the
  covariance matrix.
  (See `PR #50 <https://github.com/jrbourbeau/pyunfold/pull/50>`_)
- Adds check for negative values in ``prior`` input to ``iterative_unfold``.
  (See `PR #63 <https://github.com/jrbourbeau/pyunfold/pull/63>`_)

**Changes**:

- Makes ``data``, ``response``, ``efficiencies``, and associated uncertainties
  parameters in ``iterative_unfold`` keyword arguments. Now they can be input
  in any order. (See `PR #46 <https://github.com/jrbourbeau/pyunfold/pull/46>`_)
- Changes default prior in ``iterative_unfold`` from Jeffreys prior to a
  uniform prior. (See `PR #52 <https://github.com/jrbourbeau/pyunfold/pull/52>`_)

**Bug Fixes**:

- Fixes bug in calculation of Jeffreys prior.
  (See `PR #52 <https://github.com/jrbourbeau/pyunfold/pull/52>`_)


Version 0.2.2 (2018-05-04)
--------------------------

**Changes**:

- Adds checks for optional PyTables testing dependency. (See `PR #43 <https://github.com/jrbourbeau/pyunfold/pull/43>`_)


Version 0.2.1 (2018-05-04)
--------------------------

**Changes**:

- Adds ``MANIFEST.in`` file to package additional files with source distribution.


Version 0.2 (2018-05-04)
------------------------

**New Features**:

- Added ``SplineRegularizer`` callback to perform smoothing at each unfolding iteration. (See `PR #27 <https://github.com/jrbourbeau/pyunfold/pull/27>`_)
- Added ``Logger`` callback to write test statistic information to stdout. (See `PR #18 <https://github.com/jrbourbeau/pyunfold/pull/18>`_)
- Added contributing guide to the documentation page. (See `PR #31 <https://github.com/jrbourbeau/pyunfold/pull/31>`_)
- Added ``setup_prior`` convenience function to help validate user prior input (See `PR #35 <https://github.com/jrbourbeau/pyunfold/pull/35>`_)

**Changes**:

- Removed ``datadist.py`` module. (See `PR #19 <https://github.com/jrbourbeau/pyunfold/pull/19>`_)
- Vectorized ``covmatrix.py`` module. (See `PR #24 <https://github.com/jrbourbeau/pyunfold/pull/24>`_)
- Removed ``covmatrix`` and ``loadstats`` modules. (See `PR #38 <https://github.com/jrbourbeau/pyunfold/pull/38>`_)

**Bug Fixes**:

- Fixed typos in the online documentation (See `PR #26 <https://github.com/jrbourbeau/pyunfold/pull/26>`_)


Version 0.1 (2018-04-05)
------------------------

**New Features**:

- Initial release that includes ``iterative_unfold`` function.
