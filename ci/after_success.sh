#!/bin/bash

if [[ $COVERAGE == 'true' ]]; then
    pip install codecov
    codecov
fi
