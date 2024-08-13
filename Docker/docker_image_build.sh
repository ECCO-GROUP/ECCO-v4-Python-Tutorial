#!/bin/bash

# build docker image with build context in parent user directory
# and Dockerfile in ECCO-v4-Python-Tutorial/Docker subdirectory,
# passing current user info as build arguments


cd /home/${USER}

docker build . \
    --build-arg NB_USER=${USER} \
    --build-arg NB_UID=$(id -u ${USER}) \
    -t ecco_tut_image \
    -f ./ECCO-v4-Python-Tutorial/Docker/Dockerfile
