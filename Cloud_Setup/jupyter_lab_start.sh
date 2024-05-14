#!/bin/bash

# Shell script to start a Jupyter lab session
# in a tmux window (so it persists even if ssh tunnel is disconnected)

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'


conda init
source ~/miniforge3/bin/activate
conda activate jupyter

echo "Starting Jupyter lab session!"
echo ""

read -p 'Which port do you want to use [default 9889]: ' user_port
user_port=${user_port:-9889}
echo "using port ${user_port}"

echo ""

# Start configuration for Jupyter lab
echo "Enter password to access Jupyter lab from browser,"
echo "or leave blank to not require a password."

PW="$(python3 -c 'from jupyter_server.auth import passwd; import getpass; print(passwd(getpass.getpass(), algorithm="sha256"))')"
jlab_start="jupyter Space lab Space --no-browser Space --autoreload Space --port=\"${user_port}\" Space --ip='127.0.0.1' Space --NotebookApp.token='' Space --NotebookApp.password=\"$PW\" Space --notebook-dir=\"~/ECCO-v4-Python-Tutorial\""

# Start new tmux session
tmux new -d -s jupyterlab

# Execute commands in tmux window using send-keys
tmux send-keys -t jupyterlab source Space ~/conda/bin/activate Enter
tmux send-keys -t jupyterlab conda Space activate Space jupyter Enter
tmux send-keys -t jupyterlab ${jlab_start} Enter

# Print info about tmux session
echo -e "${red_start}Started Jupyter lab in tmux session jupyterlab"
echo -e "${red_start}Access from your local machine in a browser window at"
echo -e "${blue_start}http://127.0.0.1:${user_port}/ ${red_start}or ${blue_start}http://localhost:${user_port}/"
echo -e "${red_start}tmux session can be accessed with"
echo -e "${blue_start}tmux a -t jupyterlab"
echo -e "${red_start}and detached from current window by pressing keys"
echo -e "${blue_start}Ctrl-b d"
echo -e "${red_start}and terminated with"
echo -e "${blue_start}tmux kill-ses -t jupyterlab${nocolor_start}"
