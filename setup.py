#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

NAME = 'pyunfold'
DESCRIPTION = 'Python package for iterative Bayesian unfolding'
MAINTAINER = 'James Bourbeau'
MAINTAINER_EMAIL = 'james@jamesbourbeau.com'
URL = 'https://github.com/jrbourbeau/pyunfold'
LICENSE = 'MIT'

here = os.path.abspath(os.path.dirname(__file__))

# Want to read in package version number from __version__.py
about = {}
with open(os.path.join(here, 'pyunfold', '__version__.py'), 'r') as f:
    exec(f.read(), about)
    VERSION = about['__version__']

with open(os.path.join(here, 'requirements', 'default.txt')) as f:
    INSTALL_REQUIRES = [l.strip() for l in f.readlines() if l]

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    url=URL,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    license=LICENSE,
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
)
