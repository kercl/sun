# $\mathfrak{su}(n)$ Toolkit

This software package aims to provide several tools for working with representations of the Lie-algebras $\mathfrak{su}(n)$.

The core library is written in C for efficiency and has bindings for numpy, sympy and Mathematica.

## Setup

Python requirements:

* numpy
* sympy

## Planned Usage

### Python

Irreducibe matrix representations based on a given Dynkin labels can easily be generated. As a sample let's construct Pauli's matrices as an irrep of <sub><sub><img alt="\mathfrak{su}(2)" src="images/fa4efd26c491dfa3da955fb46c8ac023.svg" /></sub></sub> and Gell-Mann matrices as an irrep of <sub><sub><img alt="\mathfrak{su}(3)" src="images/d5f0a6fccf11f2b2a389ee9011ddd658.svg" /></sub></sub>.

```python
import sun.numeric as sun

irrep = sun.Irrep(dynkin=[1])
irrep
```
Output: 

<p align="center"><img alt="D(1)" src="images/ad0777e189b4b22f46807d3878e4e72c.svg" /></p>



This does not yet generate the matrices. However once we iterate over `irrep` or access an element throught `__getitem__` the matrices are built and stored within the object. The construction can also manually invoked through the function `construct_matrices`.

```python
irrep.construct_matrices()
for X in irrep:
    print(2 * X, end=" ")
```
Output:

<p align="center"><img alt="\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},\ \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix},\ \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}" src="images/e89d0222b50b5ffd1c309b51bc20be51.svg" /></p>

### Mathematica

```mathematica
Needs["SUN`"]

irrep = Irrep[1,0]

X = SUNLieMatrices[irrep, Method→RaisingBasis]
```

## Planned features

- [x] Generate irreducible representations from Dynkin labels
- [x] Python bindings
- [x] Mathematica bindings
- [ ] Decompose tensor products into irreducible representations


