# Dockerfile to build SMCE OSS environment for running ECCO tutorials.
FROM public.ecr.aws/smce/smce-images:smce-oss-earth-r-julia-48aae375

ENV NB_USER=jovyan
ENV HOME=/home/${NB_USER}

# USER ${NB_USER}
WORKDIR ${HOME}


# retrieve ECCO Python tutorials
RUN git clone https://github.com/ECCO-GROUP/ECCO-v4-Python-Tutorial.git

# install any additional Python packages needed
RUN mamba install -y -n notebook matplotlib ecco_v4_py

# install packages needed to run Julia notebooks
RUN julia -e 'import Pkg; Pkg.add("Pluto"); using Pluto'
RUN julia -e 'import Pkg; Pkg.add("OceanStateEstimation")'

# create symlink from home directory to jupyter_lab_start_docker.sh
RUN ln -s ${HOME}/ECCO-v4-Python-Tutorial/Docker/smce/jupyter_lab_start_smce.sh \
            ${HOME}/jupyter_lab_start_smce.sh

EXPOSE 8888

# start jupyter lab inside the container
ENTRYPOINT ["./jupyter_lab_start_smce.sh"]
