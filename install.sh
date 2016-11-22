#!/bin/bash

# Based on procedure given in https://www.scipy.org/install.html

echo "Installing..."
sudo apt-get install python-pip
sudo python -m pip install --upgrade pip
sudo pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
