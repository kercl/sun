# <sub><img alt="SU(n)" src="images/a9aa0eebca08632942699db706325eb8.svg" /></sub> Toolkit

This software package aims to provide several tools for working with representations of the spection unitary group <sub><sub><img alt="SU(n)" src="images/d0203de77fa58f64279192e742af1548.svg" /></sub></sub> and its Lie-algreba <sub><sub><img alt="\mathfrak{su}(n)" src="images/c6cc9a8850361dbf2f4bf51765b8a92d.svg" /></sub></sub>.

The core library is written in C for efficiency and has python binding for both numpy and sympy as well as Mathematica.

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
Output: <sub><sub><img alt="D(1)" src="images/4d048db56d0e2e4ff4be1777263c483e.svg" /></sub></sub>

This does not yet generate the matrices. However once we iterate over `irrep` or access an element throught `__getitem__` the matrices are built and stored within the object. The construction can also manually invoked through the function `construct_matrices`.

```python
irrep.construct_matrices()
for X in irrep:
    print(X, end=" ")
```
Output: <sub><sub><img alt="\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},\ \begin{pmatrix} 0 & i \\ -i & 0 \end{pmatrix},\ \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}" src="images/5bcaddbe725ea06ecf87ac71d47a2e08.svg" /></sub></sub>

### Mathematica

## Planned features

- [x] Generate irreducible representations from Dynkin labels
- [x] Python bindings
- [ ] Mathematica bindings
- [ ] Decompose tensor products into irreducible representations
- [ ] Matrix exponentials
- [ ] Pedagogical tools for visualizations