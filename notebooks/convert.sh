#!/bin/bash

jupyter nbconvert *ipynb --to rst
jupyter nbconvert *ipynb --to latex
jupyter nbconvert *ipynb --to html --template full
pdflatex *tex

