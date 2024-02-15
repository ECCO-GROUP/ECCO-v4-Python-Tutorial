#!/bin/bash

# # Script to enable ssh connections,
# # on EC2 instances where they are disabled by default
# # (e.g., with JPL AMIs).
# # 
# # Once this script runs successfully,
# # it should be possible to login to the instance using ssh, e.g.:
#
# ssh -i "~/.ssh/your_key_pair.pem" jpluser@your_private_ip_address


# Become root
sudo -s
 
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
        if [ $? -eq 0] ; then
            echo "Enabled sshd successfully"
        else
            echo "Error: symlink deletion did not allow sshd to be enabled"
            exit 1
    else
        echo "Error: sshd not enabled successfully"
fi
 
# Create symlink to the service (if it does not already exist)
ln -s /etc/systemd/system/multi-user.target.wants/sshd.service /usr/lib/systemd/system/sshd.service
echo "Created symlink to sshd.service"
 
# create new ssh keys
ssh-keygen -q -N "" -t rsa -b 4096 -f /etc/ssh/ssh_host_rsa_key
echo "Created new ssh keys"
 
# start sshd service
systemctl start sshd
