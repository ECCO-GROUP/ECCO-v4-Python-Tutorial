#!/bin/bash

# script to install wget and tmux 
# wget is needed to download Miniforge
# tmux is useful (but not really necessary) for running jupyter 

# # Start body of script

red_start='\033[0;31m'
blue_start='\033[0;34m'
nocolor_start='\033[0m'

if command -v dnf &> /dev/null
then
  # install wget and tmux
  sudo dnf install wget tmux -y
elif command -v apt &> /dev/null; then
  sudo apt install wget tmux -y
elif command -v yum &> /dev/null; then
  sudo yum install wget tmux -y
fi

if command -v wget &> /dev/null; then
  echo -e "${red_start}Installed wget${nocolor_start}"
else
  echo -e "${blue_start}Could not install wget! -- install it yourself ${nocolor_start}"
fi

if command -v tmux &> /dev/null; then
  echo -e "${red_start}Installed tmux${nocolor_start}"
else
  echo -e "${blue_start}Could not install tmux! -- install it yourself ${nocolor_start}"
fi

