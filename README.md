# Kronecker's Canonical Form
SageMath implementation of an algorithm to calculate Kronecker's canonical form over an exact ring.

# Setup instructions
Setup Python and SageMath on your machine, you can test it is working
by running
```sh
$ sage
```
on your terminal and see an output similar to
```sh
$ sage
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version x.y, Release Date: yyyy-mm-dd                     │
│ Using Python 3.10.z. Type "help()" for help.                       │
└────────────────────────────────────────────────────────────────────┘
sage:
```

Setup Python to handle SageMath code.
```sh
pip install sagemath
```
Test it is properly working: this is an example of a session where it is,
your output should be similar.
```sh
$ python
Python 3.10.6 (main, Aug  2 2022, 00:00:00) [GCC 12.1.1 20220507 (Red Hat 12.1.1-1)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from sage.all import *
>>> 
```

Install pytest.
```
pip install pytest
```

Install make.
# How to test
From the root of this project, run
```
make
```