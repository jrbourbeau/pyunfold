# PyUnfold

[![Build Status](https://travis-ci.org/jrbourbeau/pyunfold.svg?branch=master)](https://travis-ci.org/jrbourbeau/pyunfold)
[![codecov](https://codecov.io/gh/jrbourbeau/pyunfold/branch/master/graph/badge.svg)](https://codecov.io/gh/jrbourbeau/pyunfold)
![pypi version](https://img.shields.io/pypi/v/pyunfold.svg 'pypi version')
![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)
![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)
![license](https://img.shields.io/pypi/l/pyunfold.svg 'license')


PyUnfold is a Python implementation the D’Agostini iterative unfolding method outlined in

> G. D'Agostini, “A Multidimensional unfolding method based on Bayes' theorem”, Nucl. Instrum. Meth. A **362** (1995) 487


## Installation

The latest development version of PyUnfold can be installed directly from GitHub

```
pip install git+https://github.com/jrbourbeau/pyunfold.git
```

or you can [fork](https://guides.github.com/activities/forking/) the [GitHub repository](https://github.com/jrbourbeau/pyunfold) and install PyUnfold on your local machine via

```
git clone https://github.com/<your-username>/pyunfold.git
pip install pyunfold
```

## Quickstart

```python
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
unfolded = iterative_unfold(data, data_err,
                            response, response_err,
                            efficiencies, efficiencies_err)
```
The returned unfolded result is a dictionary containing:
- `unfolded`: Unfolded cause distribution
- `stat_err`: Statistical (Poisson) errors
- `sys_err`: Systematic errors associated with limited statistics in the response matrix

```python
print(unfolded)             
{'unfolded': array([ 94.48002622, 155.51997378]),
 'sys_err': array([0.66204237, 0.6620424 ]),
 'stat_err': array([11.2351567 , 13.75617997])}
```

## License

[MIT License](LICENSE)

Copyright (c) 2018 James Bourbeau, Zigfried Hampel-Arias
