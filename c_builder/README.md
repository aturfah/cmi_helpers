# c_builder
Python-based interface to build file catalogs.

## To install Python 2.7
- Most Mac/Linux machines will have python already installed. Check the version of python on your machine by running `python --version` from the terminal.
1. Go to [python.org](python.org) and download the Python 2.7 installer.
2. Run the installer
    - If on a Windows machine, make sure to check the "Add Python to PATH" box, this will save you a lot of grief later.
3. From the terminal, if you type `python --version` and do not get an error, then it worked.

## To install pip
- Prerequisites: Python 2.7
1. Before proceeding, from the terminal run `pip --version`. If you do not get an error, you have pip installed, and can skip this step.
2. Download `get-pip.py` from [here](https://bootstrap.pypa.io/get-pip.py)
3. Using the terminal, navigate to the directory with the `get-pip.py` file and run `python get-pip.py`
    - For macs, if you get a Permission Error, try `sudo python get-pip.py`
    - For windows users, try opening your command prompt/powershell window as Administrator, and reissue the `python get-pip.py` command.
4. Reference website [here](https://pip.pypa.io/en/stable/installing/)

## To install necessary files to run the program
- Prerequisites: Python 2.7, pip
1. Navigate to the directory with the `c_builder` directory in it.
2. Run `pip install -r requirements.txt` from the terminal.
    - If you get a Permission Error, try `pip install -r requirements.txt --user`

## To run the program
1. From the terminal, run the `run_server` script with `./run_server`
    - If on Windows use `run_server.ps1`, which can be run with `run_server.ps1`
2. Navigate to `localhost:5000` in your web browser.
