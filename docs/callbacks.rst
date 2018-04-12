.. callbacks:

:github_url: https://github.com/jrbourbeau/pyunfold

*********
Callbacks
*********

Callbacks are objects that perform some action during each iteration of
unfolding. A list of Callbacks can be passed to the
``pyunfold.iterative_unfold`` function via the ``callbacks`` keyword argument.

.. code-block:: python

    import pyunfold

    # Setup Callbacks
    logger = pyunfold.Logger()
    regularizer = pyunfold.SplineRegularizer(smooth=0.75)

    # Perform iterative unfolding
    unfolded = pyunfold.iterative_unfold(data, data_err,
                                         response, response_err,
                                         efficiencies, efficiencies_err,
                                         callbacks=[logger, regularizer])

Logger
------

Writes test statistic information for each unfolding iteration to ``sys.stdout``.
Example output:

.. code-block:: shell

    Iteration 1: ts = 0.0910, ts_stopping = 0.01
    Iteration 2: ts = 0.0476, ts_stopping = 0.01
    Iteration 3: ts = 0.0206, ts_stopping = 0.01
    Iteration 4: ts = 0.0083, ts_stopping = 0.01

.. autoclass:: pyunfold.callbacks.Logger
    :no-undoc-members:


SplineRegularizer
-----------------

.. autoclass:: pyunfold.callbacks.SplineRegularizer
    :no-undoc-members:
