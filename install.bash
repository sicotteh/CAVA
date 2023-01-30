#!/bin/bash
# Prior to running the installation script, initialize the virtual environment
# using the following style commands.(If you put the commands in a script, you will
#  have to source the script rather than execute it)
# cd venv
# export LD_LIBRARY_PATH="path_to_local_python/lib:${LD_LIBRARY_PATH}"
# source bin/activate

set -x
pip uninstall pycurl
export PYCURL_SSL_LIBRARY=nss
pip install --compile --install-option="--with-nss" --no-cache-dir pycurl
python setup.py install --user
