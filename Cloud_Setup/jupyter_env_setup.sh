#!/bin/bash

# install wget and tmux 
sh ./install_wget_tmux.sh

# download and install miniforge/conda/mamba 
sh ./install_conda_mamba.sh

# download and install a bunch of packages
# including ecco_v4_py and jupyterlab
sh ./install_packages_bigmem.sh
