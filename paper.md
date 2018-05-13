---
title: 'PyUnfold: A Python package for iterative unfolding'
tags:
  - Python
  - statistics
  - unfolding
  - deconvolution
authors:
  - name: James Bourbeau
    orcid: 
    affiliation: "1" 
  - name: Zigfried Hampel-Arias
    orcid: 0000-0003-0253-9117
    affiliation: "1, 2"
affiliations:
 - name: University of Wisconsin - Madison, Wisconsin, USA
   index: 1
 - name: IIHE, Université Libre de Bruxelles, Bruxelles, Belgium
   index: 2
date: 05 May 2018
bibliography: paper.bib
---

# Summary

``PyUnfold`` is an extensible framework for the deconvolution of discrete probability 
distributions via the iterative unfolding method described in [@agostini]. 
This method falls into the general class of inverse problems, and is especially powerful 
when no analytic form of the inverse function is explicitly available, instead requiring
an estimate (e.g. via a finite amount of simulation) of the response function.
Including the fact that the measured data comprise a finite sample, this unfolding package 
also implements the uncertainty contributions stated in [@adye2].
Full mathematical derivations of the statistical and systematic uncertainty propagation 
can be found in the accompanying LaTeX documentation.


The primary purpose of ``PyUnfold`` is to provide a toolkit for unfolding beyond the 
high-energy physics (HEP) community, whose deconvolution packages maintain a strong 
dependence on the ``ROOT`` data analysis framework [@root], almost exclusively used in HEP.
Instead, ``PyUnfold`` accepts all input objects (measured distribution, response matrix,
and associated uncertainties) as ``Numpy`` arrays, thus extending the method''s scope
to a general scientific audience.


Indeed, the method itself is data-agnostic, referring to the measurement process
as the smearing of a set of true causes into a set of detectable effects.
For example one could define as causes the true energy of a particle and the effects
the reconstructed energy of that particle in a detector.
Another example might be a set of diseases (causes) and possible clinical symptoms (effects).
So long as it is possible to encode estimable resolutions and biases connecting causes to 
effects in a binned response matrix, one can perform a deconvolution with ``PyUnfold``. 


``PyUnfold`` has been designed to be both easy to use for first-time users as well as 
flexible enough for fine-tuning an analysis.
Its main strength is its simplicity, requiring a single line to run via the ``iterative_unfold`` 
method, taking as input the user provided distributions and other runtime parameters.
This also facilitates examining the effects of varying optional parameters on the final 
unfolded distribution, allowing users to seamlessly test the robustness of results.


One such option is the initial prior probability distribution, which by default is the 
uniform prior that assigns equal probability to all possible causes.
Additionally, the user may provide a custom prior distribution as a ``Numpy`` array.
The authors have also included a utility function to calculate the non-informative 
Jeffreys'' prior.


The stopping criteria is based on a test statistic calculated by comparing unfolded 
distributions from one iteration to the next.
The user can define the strength of this criterion, as well as choose from the 
following test statistic calculations: Kolmogorov-Smirnov[@kolmogorov, @smirnov],
$\Chi^2$, relative difference, Bayes factor [@pfendner]).


Univariate spline regularization is also implemented as a means of ensuring that unfolded 
distributions do not suffer from growing fluctuations potentially arising from the finite 
binning of the response matrix.
The degree and strength of the spline can be changed from the defaults by the user.
The user can also define independent subsets of causes, which are regularized in their 
respective groups or blocks, making possible the generalization to multivariate unfolding.
Special care must be taken when performing a block unfolding, which is demonstrated in 
one of the usage examples.


Finally, for those with special data needs, the form of the covariance matrices for both
the data and response uncertainties can be changed from the defaults, with choices
being either Poisson or multinomial.
Further mathematical details can be found in the LaTeX documentation.


# Acknowledgements

The authors acknowledge support from the Wisconsin IceCube Particle Astrophysics Center
at the UW-Madison, and especially for the guidance of Stefan Westerhoff.
We also ackowledge the financial support provided by the Belgian American Educational 
Foundational Fellowship, and Wallonie-Bruxelles International.

# References
