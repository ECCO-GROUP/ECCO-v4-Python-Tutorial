# syntax=docker/dockerfile:1
FROM redhat/ubi8:latest

# clone repository and set up in container
RUN dnf update && dnf install git -y
RUN cd ~ && git clone https://github.com/ECCO-GROUP/ECCO-v4-Python-Tutorial.git
COPY ~/ECCO-v4-Python-Tutorial /

# expose port to be used later with Jupyterlab
EXPOSE 9889

# run jupyter set up script
CMD ["chmod", "755", "/ECCO-v4-Python-Tutorial/Cloud_Setup/jupyter_env_setup.sh"]
CMD ["/ECCO-v4-Python-Tutorial/Cloud_Setup/jupyter_env_setup.sh"]

## note: might have to work out how user can input Earthdata credentials when jupyter_env_setup.sh is run, and port number/password for jupyter_lab_start.sh

# run jupyter lab start script
CMD ["~/jupyter_lab_start.sh"]
