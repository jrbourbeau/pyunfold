all: install

install:
	python setup.py install

test:
	python -m pytest pyunfold

test-coverage:
	python -m pytest --cov=pyunfold --cov-report term pyunfold

test-coverage-html:
	python -m pytest --cov=pyunfold --cov-report html pyunfold
