# An Interactive Introduction to the Finite Element Method


[Joachim Schöberl](https://www.asc.tuwien.ac.at/~schoeberl)

TU Wien, [Institute of Analysis and Scientific Computing](https://www.tuwien.at/en/mg/asc)

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
element package [Netgen/NGSolve](https://www.ngsolve.org), which can be conveniently used via its
Python frontend in jupyter-notebooks.

This lecture is given in this form the first time in summer term 24.
If you have suggestions for improvements, or found some errors, please send them per mail
to the author.
Many section are still in draft version, and will be cleaned as the class proceeds.

If you like the material, show it by giving a star on github.


## Literature

* Lecture notes [Schöberl](numpde.pdf) and Faustmann+Schöberl (available in TU-Wien TUWEL)

Books:

* D.Braess: Finite Elements. Theory, Fast Solvers, and Applications in Solid Mechanics
* C. Johnson: Numerical solution of partial differential equations by the finite element method
* D.Boffi, F.Brezzi, M.Fortin: Mixed Finite Element Methods and Applications
* S.Brenner, R.Scott: The Mathematical Theory of Finite Element Methods
* A. Ern, J.-L.Guermond: Finite Elements I-III


## Installing NGSolve

Install a recent Python. Then it should be easy to install NGSolve using

    pip install jupyter numpy scipy matplotlib
    pip install --pre ngsolve
    pip install webgui_jupyter_widgets


To check the installation of NGSolve run in the console:

    python3 -c "import ngsolve; print(ngsolve.__version__)"

Then, open jupyter-notebook (or jupyter-lab or VS Code), create a new notebook, create and execute a cell with

    from ngsolve import *
    from ngsolve.webgui import Draw
    Draw (unit_cube.shape);


Known issues are
- Use pip3 instead of pip if there is no pip
- If you get an error like `externally-managed-environment`, then either use
virtual environments, or add the flag `--break-system-packages` to the pip command, see [explanation](https://veronneau.org/python-311-pip-and-breaking-system-packages.html)

- If you have conflicts with other packages, you may install NGSolve in a [virtual environment](https://docs.python.org/3/library/venv.html#creating-virtual-environments). For example I did

      python3 -m venv /Users/joachim/numpde
      source /Users/joachim/numpde/bin/activate

- If NGSolve compuatations are working, but you don't get the rendering: For jupyter notebook version < 7.0.0 you have to run additionally

      jupyter nbextension install --user --py webgui_jupyter_widgets
      jupyter nbextension enable --user --py webgui_jupyter_widgets
  


If local installation does not work, there are alternatives:

- login to a jupyter server from your browser:

  [jupyterhub.cerbsim.com](https://jupyterhub.cerbsim.com) <br>
  user: **ngshub_xx** <br>
  pwd:  **solve!xx** <br>
  with xx number from 01 to 31

  

- run NGSolve online within jupyter-lite:

  [![lite-badge](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://jschoeberl.github.io/iFEM-lite/lab?path=iFEM.ipynb)

  [https://jschoeberl.github.io/iFEM-lite/lab?path=iFEM.ipynb](https://jschoeberl.github.io/iFEM-lite/lab?path=iFEM.ipynb)

  The first time it might take a few minutes to start, and then again to import ngsolve.
  



```{tableofcontents}
```
