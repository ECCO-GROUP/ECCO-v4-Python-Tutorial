#!/bin/bash

# get path of directory where script is located
script_dir=$(dirname "$0")

# install wget and tmux 
sh ${script_dir}/install_wget_tmux.sh

# download and install miniforge/conda/mamba 
sh ${script_dir}/install_conda_mamba.sh

# download and install a bunch of packages
# including ecco_v4_py and jupyterlab
sh ${script_dir}/install_packages_bigmem.sh
