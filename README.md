## asymmetric_uncertainties



```asymmetric_uncertainties``` is a package for handling measured values with asymmetric uncertainties given by their confidence level in which the probability distribution is unknown.

Several python packages for the propagation of uncertainties are available out there,  [uncertanties package](https://pythonhosted.org/uncertainties/) being the most famous of then, however, these packages usually assume either a local linear expansion - which is only true for small uncertainties, rarely the case in astronomy - and/or assume the errors to be gaussian-distributed - therefore symmetric in respect to the nominal value, which usually is not true in high-energy astrophysics. Another characteristic of these codes is that the errors are usually the 1![equation](https://latex.codecogs.com/gif.latex?%5Csigma) values (68%) while in X-rays astronomy the errors are usually reported as the 90% confidence level value. This package then was created to solve these problems, making easy - under some assumptions - to propagate of these asymmetric uncertainties.

The values and their respective uncertainties are sampled with a pre-defined likelihood function, called "Linear Variable Width Gaussian", 
proposed in [R. Barlow's 2004 paper "Asymmetric Statistical Errors"](https://arxiv.org/abs/physics/0406120).

The operations are performed by Monte Carlo simulations, simple operations (+, -, *, /, **) are straightforward and the package has support for generic functions to be created by the user.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install ```asymmetric_uncertainties```.

```bash
pip install git+https://github.com/muryelgp/asymmetric_uncertainties.git
```

## How to use the code

```bash
from asymmetric_uncertainties import aufloat
```
A measurment with a nominal value of 10 and errors +1 and -0.5 at 68% confidence level or: 

![equation](https://latex.codecogs.com/gif.latex?a&space;=&space;10^{&plus;1.0}_{-0.5})

can be instanciated by:

```bash
a=aufloat(10, 0.5, 1, confidence=68)
```
in the same way:

![equation](https://latex.codecogs.com/gif.latex?b&space;=&space;30^{&plus;3.0}_{-3.5})


```bash
b=asymed(30, 3.5, 3, confidence=68)
```
Simple operations are straightforward, for example: 

```bash
c = a + b 
```

Check out this [jupyter notebook](https://github.com/muryelgp/PCI/blob/master/pci/How_to.ipynb) for details and more interesting examples.
