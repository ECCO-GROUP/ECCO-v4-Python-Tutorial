#!/bin/bash
jupyter nbconvert ../Tutorials_as_Jupyter_Notebooks/*ipynb --to python 
mv ../Tutorials_as_Jupyter_Notebooks/*.py .
