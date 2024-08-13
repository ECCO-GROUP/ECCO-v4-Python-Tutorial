#!/bin/bash

# Shell script to start a Jupyter lab session
# in a tmux window (so it persists even if ssh tunnel is disconnected)

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'


# Add Earthdata credentials to ~/.netrc if they are not already present
sh ./ECCO-v4-Python-Tutorial/Docker/earthdata_auth_docker.sh

echo "Starting Jupyter lab session!"
echo ""

read -p 'Which container port do you want to use [default 8888]: ' user_port
user_port=${user_port:-8888}
echo "using port ${user_port}"
echo ""

# Start configuration for Jupyter lab
echo "Enter password to access Jupyter lab from browser,"
echo "or leave blank to not require a password."

PW="$(python3 -c 'from jupyter_server.auth import passwd; import getpass; print(passwd(getpass.getpass(), algorithm="sha256"))')"


# Start Jupyter lab
source /srv/conda/bin/activate
conda activate jupyter
jupyter lab --no-browser --autoreload --port=${user_port} --ip='0.0.0.0' --NotebookApp.token='' --NotebookApp.password=\"$PW\" --notebook-dir=\"./ECCO-v4-Python-Tutorial\"

# Print info about session
echo -e "${red_start}Started Jupyter lab in Docker container"
echo -e "${red_start}Access from your local machine in a browser window"
echo -e "${red_start}by tunneling to host machine at the port"
echo -e "${red_start}you are using to connect to the Docker container, e.g., "
echo -e "${blue_start}ssh -i ~/.ssh/ec2_auth.pem -L 8888:localhost:8888 ec2-user@100.104.70.127"
echo -e "${red_start}and opening the URL in your browser"
echo -e "${blue_start}http://127.0.0.1:8888/ ${red_start}or ${blue_start}http://localhost:8888/"
echo -e "${red_start}If you are already running Jupyter on your local machine"
echo -e "${red_start}on port 8888, you need to use a different local port, e.g.,"
echo -e "${blue_start}ssh -i ~/.ssh/ec2_auth.pem -L 9889:localhost:8888 ec2-user@100.104.70.127${nocolor_start}"
