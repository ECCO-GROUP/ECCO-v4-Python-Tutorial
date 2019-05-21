*****************
Tutorial Overview
*****************


What is the format of the tutorials?
------------------------------------

All tutorials are provided as Python codes in two formats: (1) standard Python (\*.py files) and (2) Python Jupyter Notebooks (\*.ipynb files).

**Tutorials as Jupyter Notebooks**: https://github.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/tree/master/Tutorials_as_Jupyter_Notebooks

**Tutorials as Python Files**: https://github.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/tree/master/Tutorials_as_Python_Files

Because all of the tutorials are actually Python code, every one can be reproduced  locally on your own machine.  You'll only need to change the tutorials so that they point to the location of the state estimate output files on your disk and the location of your copy of the *ecco_v4_py* Python package (only if you did not install it with the pip Python package manager).


What are Jupyter notebooks?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the `Jupyter`_ website:

*"The Jupyter Notebook is an open-source web application that allows you to create and share documents that contain live code, equations, visualizations and narrative text. Uses include: data cleaning and transformation, numerical simulation, statistical modeling, data visualization, machine learning, and much more."*

.. _Jupyter : http://jupyter.org/

Jupyter notebooks allow Python code to be run in an interactive environment within a browser window.  Notebooks are divided in *cells* which can include text or Python code.  They are similar to the "command windows" of interactive desktop environments (e.g, Matlab) except that both the text, code, and the output of commands (textual output or figures) are kept together in the same document.  One of the best features of Jupyter notebooks is that they store all commands and output.  Loading someone else's notebook allows you to reproduce and build upon their work.  We hope you take advantage of this feature of the ECCO v4 tutorial!

If you want to see some examples of Juypter Notebooks before continuing further, here are some examples (1) `numerical intergration with the Trapezoid Rule`_, (2) `numerical integration with the Crank Nicolson method`_, (3) `symbolic calculations`_.  

.. _numerical intergration with the Trapezoid Rule: http://nbviewer.jupyter.org/github/ipython/ipython/blob/4.0.x/examples/IPython%20Kernel/Trapezoid%20Rule.ipynb
.. _symbolic calculations : http://nbviewer.jupyter.org/github/ipython/ipython/blob/4.0.x/examples/IPython%20Kernel/SymPy.ipynb
.. _numerical integration with the Crank Nicolson method : http://nbviewer.jupyter.org/github/waltherg/notebooks/blob/master/2013-12-03-Crank_Nicolson.ipynb

This `notebook basics`_ page may be helpful for learning how to navigate around in notebooks.

.. _notebook basics : http://nbviewer.jupyter.org/github/jupyter/notebook/blob/master/docs/source/examples/Notebook/Notebook%20Basics.ipynb


Will I learn Python just from reading these tutorials?
------------------------------------------------------

Unlikely!  These tutorials are not a comprehensive introduction to Python.  **Nevertheless, throughout the tutorial there are many helpful tips to guide the non-native Python speaker.** These tips come early on so it helps to start from the beginning.  


What Python should I review before getting started?
---------------------------------------------------

There are thousands of Python packages, many of which would no doubt be useful for working with the ECCO v4 state estimate.  Howefver, there are *three* that you will should familiarze yourself with on some level before starting: **NumPy**, **Matplotib**, and **xarray**.  **NumPy** provides n-dimensional matrices and matrix operations, **Matplotlib** provides plotting tools, and **xarray** provides a framework for labelling and easily accessing subsets of the **NumPy** n-dimensional matrices.  **xarray** is not necessary for working with ECCO v4 output as one does not need to have labeled dimensions and coordinatess to analyze arrays.  Nevertheless, you may find that working with arrays that have well  labeled coordinates and dimensions makes life much, much easier.

NumPy
^^^^^
From the `NumPy website <http://www.numpy.org/>`_

NumPy is the fundamental package for scientific computing with Python. It contains among other things:
    a powerful N-dimensional array object
    sophisticated (broadcasting) functions
    tools for integrating C/C++ and Fortran code
    useful linear algebra, Fourier transform, and random number capabilities


Matplotlib
^^^^^^^^^^
From the `Matplotlib website <https://matplotlib.org/>`_

Matplotlib is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms. Matplotlib can be used in Python scripts, the Python and IPython shell, the jupyter notebook, web application servers, and four graphical user interface toolkits.

Matplotlib tries to make easy things easy and hard things possible. You can generate plots, histograms, power spectra, bar charts, errorcharts, scatterplots, etc., with just a few lines of code.

For simple plotting the pyplot module provides a MATLAB-like interface, [...]. For the power user, you have full control of line styles, font properties, axes properties, etc, via an object oriented interface or via a set of functions familiar to MATLAB users.

xarray
^^^^^^
From the `xarray website <http://xarray.pydata.org/en/stable/why-xarray.html>`_

Adding dimensions names and coordinate indexes to numpy's ndarray_ makes many powerful array operations possible.

The N-dimensional nature of xarray's data structures makes it suitable for dealing
with multi-dimensional scientific data, and its use of dimension names instead of axis labels (``dim='time'`` instead of ``axis=0``) makes such arrays much more manageable than the raw numpy ndarray: [...]

.. _ndarray : https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

What if I don't like the way you do X?
--------------------------------------

That's ok.  We're always open to suggestions.  


What if I find a mistake?
-------------------------

Definitely please tell us!  You can write Ian.Fenty at jpl.nasa.gov


What if I would like to contribute with a tutorial of my own?
-------------------------------------------------------------

That would be fantastic!  We're all in this together.  If you have a notebook that you think would be helpful, let us know and we'll do our best to integrate it.  


Bonus Tutorials
---------------

Two valuable tutorials are included at the end as **Bonus Tutorials**.  Both bonus tutorials use the xmitgcm package to load MITgcm model output. Consequently, some of the synatax of their calculations will be different than the syntax used in the main tutorial.  For the reader who has completed the main tutorial the bonus tutorials will be easy to understand. 

1) "Evaluating budgets in the ECCOv4 model run using xgcm" by Erik Tesdal 

2) "Vector calculus in ECCO: The Transport, divergence, vorticity and the Barotropic Vorticity Budget" by Maike Sonnewald

