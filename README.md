# Kronecker's Canonical Form
SageMath implementation of an algorithm to calculate Kronecker's canonical form over an exact ring.

# Setup instructions
First, it is required to clone the repository.

In order to test the code provided or include it in scripts, it is required to follow the following steps.

Setup **Python** on your machine; the version the algorithm has been implemented with
is \textbf{Python 3.10.6}.

Setup **SageMath** on your machine and make sure Python can handle SageMath code, this should be achieved
by installing the Python package _sagemath_

```sh
pip install sagemath
```

On success, the following snippet does not produce any
errors.
```py
from sage.all import *
```

Last, install _make_ and the Python package _pytest_.

# How to test
From the root of this project, run
```
make tests
```

It is also possible to get a report on the test coverage provided by the test suite.
```sh
coverage run -m pytest -v tests && coverage report -m
```