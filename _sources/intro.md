# An Interactive Introduction to the Finite Element Method


[Joachim Schöberl](https://www.asc.tuwien.ac.at/~schoeberl)

TU Wien, [Institute of Analysis and Scientific Computing](https://www.asc.tuwien.ac.at)

The finite element method is a powerful tool for computer simulation of problems in engineering and sciences.
Such problems are often described mathematically by means of partial
differential equations, which are then discretized on a mesh.

Finite element methods are in the intersection of mathematics,
engineering, and scientific computing. Thus it is natural that there
are very different courses teaching finite elements, from very
theoretical courses never touching the computer, to very applied
classes running commercial programmes without teaching the methods
behind.

In this class we follow an approach in between: We aim explaining the
mathematical theory, and giving students the possibility to try all
methods on the computer. For this we are using the open source finite
element package Netgen/NGSolve, which can be conveniently used via its
Python frontend in jupyter-notebooks.

This lecture is given in this form the first time in summer term 24.
If you have suggestions for improvements, or found some errors, please send them per mail
to the author, or open a pull request on the source repo of the book.
Many section are still in draft version, and will be cleaned as the class proceeds.


## Literature

* Lecture notes Schöberl and Faustmann+Schöberl

Books:

* D.Boffi, F.Brezzi, M.Fortin: Mixed Finite Element Methods and Applications
* D.Braess: Finite Elements. Theory, Fast Solvers, and Applications in Solid Mechanics
* S.Brenner, R.Scott: The Mathematical Theory of Finite Element Methods
* C. Johnson: Numerical solution of partial differential equations by the finite element method
* A. Ern, J.-L.Guermond: Finite Elements I-III


## Installing NGSolve

Install a recent Python. Then it should be easy to install NGSolve using

    pip install jupyter numpy scipy matplotlib
    pip install ngsolve
    pip install webgui_jupyter_widgets


To check the installation of NGSolve run in the console:

    python3 -c "import ngsolve; print(ngsolve.__version__)"

Then, open jupyter-notebook, create a new notebook, create and execute a cell with

    from ngsolve import *
    from ngsolve.webgui import Draw
    Draw (unit_cube.shape);


Known issues are
- Use pip3 instead of pip if there is no pip
- If you get an error like `externally-managed-environment`, then either use
virtual environments, or add the flag `--break-system-packages` to the pip command, see [explanation](https://veronneau.org/python-311-pip-and-breaking-system-packages.html)


If local installation does not work, there are alternatives:

- login to a jupyter server from your browser:

  [jupyterhub.cerbsim.com](https://jupyterhub.cerbsim.com) <br>
  user: **ngshub_xx** <br>
  pwd:  **solve!xx** <br>
  with xx number from 01 to 31

  

- run NGSolve online within jupyter-lite:

  [![lite-badge](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https:n.ipynb)


  [https://ngsolve.github.io/jupyterlite_ngsolve/lab?path=poisson.ipynb](https://jschoeberl.github.io/iFEM-lite/lab?path=iFEM.ipynb)

  The first time it might take a few minutes to start, and then again to import ngsolve.
  



```{tableofcontents}
```
