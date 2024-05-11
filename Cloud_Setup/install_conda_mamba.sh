#!/bin/bash

# Shell script fgr setting up Miniforge/conda, jupyter, essential Python packages on an AWS EC2 instance.
# Files will go into the ~/conda/ directory 

# # Start body of script

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'


# only proceed if both tmux and wget are installed

if ! command -v wget &> /dev/null; then
   echo "wget not installed, install first then re-run this script"
   exit 1
elif ! command -v tmux &> /dev/null; then
   echo "tmux not installed, install first then re-run this script"
   exit 1
fi


# retrieve and install miniforge
echo -e "${red_start}Starting Miniforge3 installation${nocolor_start}"

# download Miniforge 
if [ -f ~/Miniforge3.sh ]; then
   wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O ~/Miniforge3.sh
else
   echo "Miniforge already downloaded!"
fi

# install Miniforge 
# default directory is ~/miniforge3

sh ~/Miniforge3.sh -b

# add conda and mamba to .bashrc or equivalent
~/miniforge3/bin/conda init
~/miniforge3/bin/mamba init

echo -e "${red_start}Completed Miniforge3 installation${nocolor_start}"
echo -e "${red_start}Restart your shell${nocolor_start}"

mamba update -n base -c conda-forge conda

# this one is a bit of a mystery, but required
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniforge3/lib" >> ~/.bashrc
