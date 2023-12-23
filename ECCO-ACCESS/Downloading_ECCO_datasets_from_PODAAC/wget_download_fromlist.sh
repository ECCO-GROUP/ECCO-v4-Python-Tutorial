#!/bin/bash

# Andrew Delman, December 2023
# 
# This script takes PO.DAAC Opendap URLs listed in a text file
# and downloads them to the current directory, or another directory
# specified by option -P.
# For example:
#     ./wget_download_fromlist.sh -i urls_download.txt \
#        -P /ECCOv4_downloads/ -n Caribbean \
#        -u username -p password
# downloads the files from URLs listed in ./urls_download.txt,
# to the directory /ECCOv4_downloads/, and appends the
# identifier 'Caribbean' to each of the downloaded file names.
# Input options can be specified either using the -i -P -n -u -p tags
# shown above, or sequentially in that order without the tags.
# However, option and positional/sequential inputs can not be combined
# when this script is called.
# 
# Note: if NASA Earthdata user authentication is already stored in
# the user's .netrc file, then -u and -p can be omitted,
# and storing the authentication in .netrc is recommended for frequent users.
# See https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
# for a step-by-step guide to set this up.


# default arguments
download_dir="./"
append_id=""

# positional argument support
if [ "$1" != "-i" ]; then
  download_dir="$2"
  append_id="$3"
  username="$4"
  password="$5"
fi

# if input options given, assign to string variables
while getopts ":i:P:n:u:p:" option; do
  case $option in
    i) # text file specifying URLs to download
       url_listfile="$OPTARG";;
    P) # directory (with path) to download files to
       download_dir="$OPTARG";;
    n) # identifier to append to file names
       append_id="$OPTARG";;
    u) # Earthdata username
       username="$OPTARG";;
    p) # Earthdata password
       password="$OPTARG";;
  esac
done

# if -i option not supplied then assume $1 is URL file list
if [ -z ${url_listfile+x} ]; then
  url_listfile="$1"
fi

# if no input arguments supplied, return error message and exit
if [ -z ${1+x} ]; then
  echo "Error: no URL file list supplied. No files downloaded."
  exit
fi

# create download directory if it does not already exist
mkdir -p "$download_dir"


# read URLs from URL text file and download to $download_dir,
# with file names as the name of the granule plus $append_id
while IFS= read -r line; do
  no_paths=${line##*/granules/}
  after_dap=${no_paths#*.dap.}
  filename=${no_paths%.dap.nc*}
  if [ "${after_dap:0:2}" = "nc" ]; then
    filename=$filename"_"$append_id".nc"
  else
    echo "Downloaded file type uncertain; may not be NetCDF"
  fi
  if [ "${download_dir:(-1)}" = "/" ]; then
    path_filename=$download_dir$filename
  else
    path_filename=$download_dir"/"$filename
  fi
  if [ -n "$username" ] && [ -n "$password" ]; then
    wget -nv -nc -c -O $path_filename \
     --user="$username" --password="$password" $line
  else
    wget -nv -nc -c -O $path_filename $line
  fi
done < $url_listfile
