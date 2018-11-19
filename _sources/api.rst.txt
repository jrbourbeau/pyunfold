.. _api:

:github_url: https://github.com/jrbourbeau/pyunfold

************
PyUnfold API
************

Iterative unfolding
-------------------

.. autofunction:: pyunfold.iterative_unfold


Callbacks
---------

.. autoclass:: pyunfold.callbacks.Logger
    :no-undoc-members:

.. autoclass:: pyunfold.callbacks.SplineRegularizer
    :no-undoc-members:


Priors
------

.. autofunction:: pyunfold.priors.uniform_prior

.. autofunction:: pyunfold.priors.jeffreys_prior


Test Statistics
---------------

.. autofunction:: pyunfold.teststat.get_ts

.. autoclass:: pyunfold.teststat.KS
    :members:

.. autoclass:: pyunfold.teststat.Chi2
    :members:

.. autoclass:: pyunfold.teststat.RMD
    :members:

.. autoclass:: pyunfold.teststat.BF
    :members:
