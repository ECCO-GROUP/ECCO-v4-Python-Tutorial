#!/bin/bash

# Shell script for setting up conda, jupyter, essential Python packages on an AWS EC2 instance.

# before installing, wise to update the local
# machine with the latest packages. 
# on redhat linux machines, 
# $ sudo dnf update -y
# $ sudo dnf install git -y

# # Start body of script

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'


source ~/.bashrc
conda init

# create conda environment called 'jupyter'
env_name='jupyter'

conda create --name ${env_name} python=3.12 -y
echo -e "${red_start}Created ${env_name} conda environment${nocolor_start}"

# install python packages (using mamba) in jupyter environment
conda activate ${env_name} 

echo -e "${red_start}Installing Python packages in jupyter environment${nocolor_start}"

# -n syntax installs the packages into the specified environment name
mamba install -n ${env_name} requests tqdm numpy pandas xorg-libice libexpat libevent -y
mamba install -n ${env_name} nspr alsa-lib libogg libpq xorg-renderproto xorg-xf86vidmodeproto graphite2 expat -y
mamba install -n ${env_name} libgpg-error dbus libflac gettext xcb-util-wm xorg-libx11 xcb-util-image -y
mamba install -n ${env_name} xkeyboard-config  libxkbcommon fonts-conda-forge font-ttf-ubuntu gstreamer zlib -y
mamba install -n ${env_name} xorg-xextproto libpng attr mpg123 pixman libvorbis glib-tools -y
mamba install -n ${env_name} libsystemd0 xcb-util-keysyms xorg-libxrender libllvm15 -y
mamba install -n ${env_name} font-ttf-dejavu-sans-mono pcre2 font-ttf-inconsolata font-ttf-source-code-pro -y
mamba install -n ${env_name} lame nss xorg-xproto pthread-stubs xorg-libxdmcp -y
mamba install -n ${env_name} libgcrypt xorg-libsm xorg-libxext fonts-conda-ecosystem xorg-kbproto mysql-libs -y
mamba install -n ${env_name} fontconfig libjpeg-turbo xcb-util-renderutil -y
mamba install -n ${env_name} glib freetype libcap libcups libopus -y
mamba install -n ${env_name} gst-plugins-base mysql-common xcb-util -y
mamba install -n ${env_name} cairo libsndfile harfbuzz xorg-libxau -y
mamba install -n ${env_name} libglib libxcb qt-main pyqt matplotlib netcdf4 -y
mamba install -n ${env_name} h5netcdf boto3 lxml scipy goes proj pyproj cartopy notebook -y 
mamba install -n ${env_name} progressbar gsw nco pympler -y
mamba install -n ${env_name} xarray[complete] jupyterlab dask_labextension s3fs -y
mamba install -n ${env_name} pyresample -y

# install remaining packages using pip
pip install ecco_access --no-cache-dir
pip install ecco_v4_py --no-cache-dir

echo -e "${red_start}Completed Python package installations${nocolor_start}"

# Set up NASA Earthdata credential

echo -e "${red_start}Setting up NASA Earthdata authentication${nocolor_start}"
# NASA Earthdata authentication
# check if credentials are already archived in ~/.netrc, and if not then prompt the user for them
earthdata_cred_stored=0
if [ -f ~/.netrc ]; then
    if grep -q "machine urs.earthdata.nasa.gov" ~/.netrc; then
        earthdata_cred_stored=1
        echo -e "${red_start}Earthdata credentials already archived ${nocolor_start}"
    fi
fi
if [ $earthdata_cred_stored -eq 0 ]; then
    if [ -f ~/.netrc ]; then sudo chmod 600 ~/.netrc; fi
    read -p 'NASA Earthdata username: ' uservar
    read -sp 'NASA Earthdata password: ' passvar
    echo -e "machine urs.earthdata.nasa.gov\n    login ${uservar}\n    password ${passvar}\n" >> ~/.netrc
    
    echo -e "\n${red_start}NASA Earthdata authentication info archived in ~/.netrc${nocolor_start}"
fi
sudo chmod 400 ~/.netrc

# create symlink to jupyter_lab_start.sh from the user's home directory
ln -s ~/ECCO-v4-Python-Tutorial/Cloud_Setup/jupyter_lab_start.sh ~/jupyter_lab_start.sh

echo "goodbye!"
