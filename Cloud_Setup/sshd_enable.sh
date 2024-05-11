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


# default flag on whether or not to try to 
# completely reinstall openssh
reinstall_ssh=false

# Try to enable and start sshd
systemctl enable sshd --now
if [ $? -eq 0 ]; then
    echo "Enabled and started sshd successfully"
else
    # if ssh is not enabled, then see if removing a dummy softlink to /dev/null helps
    # .... useful for certain AMIs
    if [ $( readlink -f /etc/systemd/system/sshd.service) = "/dev/null" ]; then
        # Delete this /dev/null symlink 
        rm -f /etc/systemd/system/sshd.service
        echo 'Deleted symlink to /dev/null'
         
        # Re-try enabling sshd
        systemctl enable sshd --now
        if [ $? -eq 0 ]; then
            echo "Enabled sshd successfully"
        else
            echo "Error: symlink deletion did not allow sshd to be enabled"
            reinstall_ssh=true 
        fi
    else
        # the softlink hack didn't work, sshd still doesn't work
        echo "Error: sshd not enabled successfully"
        reinstall_ssh=true 
    fi
fi


if $reinstall_ssh ; then
    # delete old key
    rm /etc/ssh/ssh_host_rsa_key --force

    if ! command -v dnf &> /dev/null ; then
      echo "dnf could not be found"
      echo "this machine probably isn't using redhat"
      echo "install openssh-server using your own methods"
      exit 1
    fi

    # uninstall openssh server
    dnf remove openssh-server -y

    # reinstall service
    dnf install openssh-server -y
fi

 

# Create symlink to the service (if it does not already exist)
if [ ! -f /usr/lib/systemd/system/sshd.service ]; then
    ln -s /etc/systemd/system/multi-user.target.wants/sshd.service /usr/lib/systemd/system/sshd.service
    echo "Created symlink to sshd.service"
fi
 
# create new ssh keys
read -p 'Do you want to generate new ssh keys (Y/[N]) ' new_keys 
new_keys=${new_keys:-N}
if [ $new_keys == "Y" ] || [ $new_keys == "y" ]; then
  echo "Generating new ssh keys ...."
  ssh-keygen -q -N "" -t rsa -b 4096 -f /etc/ssh/ssh_host_rsa_key <<< y
  echo "Created new ssh keys"
fi 

# Try to enable sshd
systemctl enable sshd --now
if [ $? -eq 0 ]; then
    echo "Enabled and started sshd successfully"
else
    echo "Despite everything, we could not get sshd to work"
    exit 1
fi


echo "Now you can login to your instance using ssh, e.g.:"
echo '$ ssh -i "~/.ssh/your_key_pair.pem" jpluser@private_ip_address'


echo ""
# move git repo from the 'ssm-user' account to the 'ssh user' account.
# and change ownership (if requested)

read -p 'Do you want to move the ECCO-v4-Python-Tutorial repo from 'ssm-user' account to a different account? (Y/[N]) ' move_opt
move_opt=${move_opt:-N}
if [ $move_opt == "Y" ] || [ $move_opt == "y" ]; then
    read -p 'User name of new owner [note: for JPL use "jpluser"]: ' ssh_user
    mkdir /home/${ssh_user}/git_repos
    cd /home
    mv ./ssm-user/ECCO-v4-Python-Tutorial ./${ssh_user}/git_repos/
    echo "Moved ECCO-v4-Python-Tutorial repo to /home/${ssh_user}/git_repos/"
    chown -R ${ssh_user}:${ssh_user} ./${ssh_user}/git_repos/ECCO-v4-Python-Tutorial
    echo "Changed owner and group of git repo to ${ssh_user}"
fi
