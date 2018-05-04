.. _changelog:

:github_url: https://github.com/jrbourbeau/pyunfold

*************
Release Notes
*************

Version 0.2.dev0 (TBD)
----------------------

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
