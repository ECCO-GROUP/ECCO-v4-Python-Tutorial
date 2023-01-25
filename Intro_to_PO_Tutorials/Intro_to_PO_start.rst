*****************
Intro to PO Tutorials: Getting Started
*****************


What are the Intro to PO Tutorials?
------------------------------------

The tutorials in this series use output from the ECCO version 4 release 4 (v4r4) state estimate to illustrate foundational concepts in the physics of the ocean (physical oceanography or PO for short). These tutorials are written as Jupyter notebooks; this format allows the concepts, code, and results of running the code to be viewed together in one document.
While these notebooks can be read online, it is strongly recommended to download them, and the ECCO output needed to run them, in order to allow users to interact with the data themselves. You can even tinker with the notebooks yourself to look at different regions or perform different calculations, that’s part of the fun!


Who are these tutorials for?
----------------------------

These tutorials are geared towards graduate students currently or recently enrolled in physical oceanography or geophysical fluid dynamics (GFD) classes.  They are intended to supplement instruction with practical application of concepts in an ocean state estimate that is actively used by ocean and climate science researchers.  Brief reviews of concepts are included with each tutorial, with references to textbooks commonly used in classrooms.

That said, everyone is welcome to try out these tutorials, and even the experienced researcher may find useful code or informative tips on how to analyze the ECCO state estimate.


What do I need to get started?
------------------------------

You will need to download Jupyter Notebook in order to interact with the tutorial notebooks.  For the user new to Python, probably the easiest way to get Jupyter Notebook is by downloading `Anaconda`_, which includes Jupyter software as well as a number of common Python packages that are used in the tutorial notebooks. 

.. _Anaconda : https://www.anaconda.com/products/distribution

While the tutorials are intended to help users get more comfortable with scientific computing in Python, a little familiarity with numeric computations in Python is recommended before diving into the Intro to PO Tutorials.  The `Python and Python Packages`_ tutorial links to some helpful resources for the Python beginner.  The tutorial also tells you how to download the *ecco_v4_py* package, which you will need for some computations in the Intro to PO tutorial series.

.. _Python and Python Packages : https://ecco-v4-python-tutorial.readthedocs.io/Installing_Python_and_Python_Packages.html

It is strongly recommended to review the first few tutorials under the **Getting Started** heading before starting the Intro to PO Tutorials.  These will tell you `what the ECCOv4 state estimate is`_, how ECCOv4 output `is structured`_, how to `get the tools`_ you need to start working with ECCOv4 output in Python, and `how to download`_ ECCOv4 output files.

.. _what the ECCOv4 state estimate is : https://ecco-v4-python-tutorial.readthedocs.io/intro.html
.. _is structured : https://ecco-v4-python-tutorial.readthedocs.io/fields.html
.. _get the tools : https://ecco-v4-python-tutorial.readthedocs.io/Installing_Python_and_Python_Packages.html
.. _how to download : https://ecco-v4-python-tutorial.readthedocs.io/Downloading_ECCO_Datasets_from_PODAAC_Python.html


Which concepts are covered?
---------------------------

The following concepts have tutorials that are either completed or planned for completion by summer 2023:

**Geostrophic Dynamics**

- Geostrophic balance
- Thermal wind
- Steric height

**Ageostrophic Dynamics**

- Ekman dynamics
- Equatorial waves

**Vorticity/Quasi-Geostrophy**

- Sverdrup balance
- Rossby waves
- Vorticity budget

**Property Budgets**

- Potential density/water masses
- Potential temperature
- Salinity, salt, and freshwater
- Sea ice

**Climate Topics**

- MOC and transports
- Earth energy imbalance
- Sea level rise
- Hydrological cycle