#!/bin/bash

# Install Python virtual environment
echo "[INFO] Set up virtual Python environment using Python 3 virtualenv module."
python3 -m virtualenv py_venv

# Activate Python virtual environment
echo "[INFO] Activate virtual environment (you have to do this in every new shell!). Now we can install packages via pip locally. You can leave the virtual environment by calling 'deactivate'."
source py_venv/bin/activate

# Upgrade pip to newest version
pip install --upgrade pip

# Install needed packages
pip install numpy scipy pandas matplotlib

echo "[INFO] Keep in mind: You can activate the virtual environment with all dependencies with 'source py_venv/bin/activate'."
