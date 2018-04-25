.. _overview:

:github_url: https://github.com/jrbourbeau/pyunfold

************
Overview
************

PyUnfold is a Python language software package for the incorporation of the imperfections
of the measurement process in the data analysis roadmap.

In an ideal world, scientists would have access to the perfect detector:
an apparatus that makes no error in measuring a desired quantity.
However, real detectors have finite resolutions and can have characteristic biases that may 
not be eliminable, yet themselves can be measured or estimated.
Statistical uncertainties also propagate through the measurement process, inclusion
of which typically requires either careful book-keeping or conservative estimates.

By building a matrix that encodes a detector's smearing of the desired true quantity
into the measured observable(s), a deconvolution can be performed that attempts to back out
an estimate of the true variable.
The unfolding method implemented in PyUnfold accomplishes this deconvolution
by harnessing the power of Bayes' Theorem in an iterative procedure, providing results
based on physical expectations of the true quantity.


------------------
Who uses PyUnfold?
------------------

The potential audience for PyUnfold includes natural scientists who are in need of a 
statistical and programming framework in which to incorporate the uncertainties of 
their measurement processes to estimate their desired variables.
PyUnfold has been successfully implemented  to unfold the cosmic-ray energy spectrum using data
from the HAWC_ and IceCube_ observatories.

.. HAWC_: https://www.hawc-observatory.org/
.. IceCube_: https://icecube.wisc.edu/


-----
Goals
-----

PyUnfold is intended to provide experimental scientists

- a tool for the analysis of measurement effects on physical distributions
- a straightforward programming interface that is suitable for many experimental applications
- the ability to painlessly incorporate and propagate statistical and systematic uncertainties



-------
History
-------

The original version of PyUnfold was brought to life by Zigfried Hampel-Arias in early 2016,
and was later diligently consolidated into its present form by James Bourbeau, both proud UW Badgers.
The first public release of PyUnfold was made available in April 2018.


