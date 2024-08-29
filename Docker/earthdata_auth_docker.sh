#!/bin/bash

# Shell script for adding NASA Earthdata authentication credentials to ~/.netrc
# if they are not already in the file

# # Start body of script

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'

# Set up NASA Earthdata credential

echo "${red_start}Setting up NASA Earthdata authentication${nocolor_start}"
# NASA Earthdata authentication
# check if credentials are already archived in ~/.netrc, and if not then prompt the user for them
earthdata_cred_stored=0
if [ -f ~/.netrc ]; then
    if grep -q "machine urs.earthdata.nasa.gov" ~/.netrc; then
        earthdata_cred_stored=1
        echo "${red_start}Earthdata credentials already archived ${nocolor_start}"
    fi
fi
if [ $earthdata_cred_stored -eq 0 ]; then
    if [ -f ~/.netrc ]; then chmod 600 ~/.netrc; fi
    read -p 'NASA Earthdata username: ' uservar
    read -p 'NASA Earthdata password: ' passvar
    echo "machine urs.earthdata.nasa.gov\n    login ${uservar}\n    password ${passvar}\n" >> ~/.netrc
    
    echo "\n${red_start}NASA Earthdata authentication info archived in ~/.netrc${nocolor_start}"
fi
chmod 400 ~/.netrc
