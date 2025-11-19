#!/bin/bash

# Shell script for setting up conda, jupyter, essential Python packages on an AWS EC2 instance.
# Assumes that the ECCO-v4-Python-Tutorial Github repository has already been downloaded using:
# 
# $ sudo dnf update -y
# $ sudo dnf install git -y
# $ cd ~
# $ git clone https://github.com/ECCO-GROUP/ECCO-v4-Python-Tutorial.git

# Then run this script:
# 
# $ sudo chmod 755 ~/ECCO-v4-Python-Tutorial/Cloud_Setup/jupyter_env_setup.sh
# $ ~/ECCO-v4-Python-Tutorial/Cloud_Setup/jupyter_env_setup.sh



# # Start body of script

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'


source ~/.bashrc
mamba init

# create jupyter environment
mamba create --name jupyter python=3.11 -y
echo -e "${red_start}Created jupyter environment${nocolor_start}"

# install python packages (using mamba) in jupyter environment
mamba activate jupyter
echo -e "${red_start}Installing Python packages in jupyter environment${nocolor_start}"
mamba install requests tqdm numpy pandas -y
mamba install xorg-libice libexpat libevent -y
mamba install nspr alsa-lib libogg libpq -y
mamba install xorg-renderproto xorg-xf86vidmodeproto graphite2 expat -y
mamba install libgpg-error dbus -y
mamba install libflac gettext -y
mamba install xcb-util-wm xorg-libx11 xcb-util-image -y
mamba install xkeyboard-config -y
mamba install libxkbcommon fonts-conda-forge font-ttf-ubuntu gstreamer zlib -y
mamba install xorg-xextproto libpng attr mpg123 -y
mamba install pixman libvorbis glib-tools -y
mamba install libsystemd0 xcb-util-keysyms xorg-libxrender libllvm15 -y
mamba install font-ttf-dejavu-sans-mono pcre2 font-ttf-inconsolata font-ttf-source-code-pro -y
mamba install lame nss xorg-xproto pthread-stubs xorg-libxdmcp -y
mamba install libgcrypt xorg-libsm xorg-libxext fonts-conda-ecosystem xorg-kbproto mysql-libs -y
mamba install fontconfig libjpeg-turbo xcb-util-renderutil -y
mamba install glib -y
mamba install freetype libcap libcups libopus -y
mamba install gst-plugins-base mysql-common xcb-util -y
mamba install cairo -y
mamba install libsndfile harfbuzz xorg-libxau -y
mamba install libglib libxcb -y
mamba install qt-main -y
mamba install pyqt -y
mamba install matplotlib -y
mamba install netcdf4 -y
mamba install h5netcdf -y
mamba install boto3 lxml -y
mamba install scipy -y
mamba install geos -y
mamba install proj pyproj -y
mamba install cartopy -y
mamba install notebook -y
mamba install progressbar -y
mamba install gsw -y
mamba install nco -y
mamba install pympler -y

# install remaining packages using pip
# (mamba installs tend to get killed on t2.micro)
pip install dask
pip install "xarray[complete]"
pip install jupyterlab
pip install dask_labextension
pip install s3fs
pip install ecco_access
pip install ecco_v4_py

echo -e "${red_start}Completed Python package installations${nocolor_start}"


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
