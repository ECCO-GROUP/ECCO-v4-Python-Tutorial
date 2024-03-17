#!/bin/bash

# Script to enable ssh connections,
# on EC2 instances where they are disabled by default
# (e.g., with JPL AMIs).
# 
# This script must be run as root, otherwise an error is returned:
# $ sudo ./sshd_enable.sh
# 
# Once this script runs successfully,
# it should be possible to login to the instance using ssh, e.g.:
#
# $ ssh -i "~/.ssh/key_pair.pem" jpluser@private_ip_address


# Return error if not running as root
if [ $( whoami ) != "root" ]; then
    echo "Error: this script must be run as root"
    echo "Please re-run using sudo, e.g.:"
    echo "$ sudo ./sshd_enable.sh"
    exit 1
fi


# Try to enable sshd
systemctl enable sshd
if [ $? -eq 0 ]; then
    echo "Enabled sshd successfully"
else
    if [ $( readlink -f /etc/systemd/system/sshd.service) = "/dev/null" ]; then
        # Delete this /dev/null symlink 
        rm -f /etc/systemd/system/sshd.service
        echo 'Deleted symlink to /dev/null'
         
        # Re-try enabling sshd
        systemctl enable sshd
        if [ $? -eq 0 ]; then
            echo "Enabled sshd successfully"
        else
            echo "Error: symlink deletion did not allow sshd to be enabled"
            exit 1
        fi
    else
        echo "Error: sshd not enabled successfully"
    fi
fi
 
# Create symlink to the service (if it does not already exist)
if [ ! -f /usr/lib/systemd/system/sshd.service ]; then
    ln -s /etc/systemd/system/multi-user.target.wants/sshd.service /usr/lib/systemd/system/sshd.service
    echo "Created symlink to sshd.service"
fi
 
# create new ssh keys
ssh-keygen -q -N "" -t rsa -b 4096 -f /etc/ssh/ssh_host_rsa_key
echo "Created new ssh keys"
 
# start sshd service
systemctl start sshd
echo "Started sshd"
echo "Now you can login to your instance using ssh, e.g.:"
echo '$ ssh -i "~/.ssh/your_key_pair.pem" jpluser@private_ip_address'


# move git repo to ssh user's directory and change ownership (if requested)
read -p 'Move ECCO-v4-Python-Tutorial repo to different user? (Y/[N]) ' move_opt
if [ $move_opt == "Y" ] || [ $move_opt == "y" ]; then
    read -p 'User name of new owner [jpluser for JPL]: ' ssh_user
    cd /home
    mv ./ssm-user/ECCO-v4-Python-Tutorial ./${ssh_user}/
    echo "Moved ECCO-v4-Python-Tutorial repo to /home/${ssh_user}/"
    chown -R ${ssh_user}:${ssh_user} ./${ssh_user}/ECCO-v4-Python-Tutorial
    echo "Changed owner and group of git repo to ${ssh_user}"
fi
