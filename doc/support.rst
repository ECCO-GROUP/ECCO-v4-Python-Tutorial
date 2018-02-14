Getting Help
============

The ECCO Support Mailing List
-----------------------------

For questions or comments please contact us via: ecco-support@mit.edu


Problems installing Python libraries?
-------------------------------------

Sometimes the pip Python package manager doesn't download the latest version of libraries from pypi, particlar on Mac/osX systems, because old versions of the libraries can be stored in the local cache directory.  Deleting the pip cache directory and always using the "--no-cache-dir" option when downloading/installing packages solves this problem. 

::

    rm -fr ~/Library/Caches/pip/*
    pip install ecco_v4_py --no-cache-dir


Sometimes libraries get updated.  To update packages installed with conda use the the following command:

::

    conda update anaconda