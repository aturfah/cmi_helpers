# c_builder
Python-based interface to build file catalogs.

## To install pip
- Prerequisites: Python 2.7
1. Before proceeding, from the terminal run `pip --version`. If you do not get an error, you have pip installed, and can skip this step.
2. Download `get-pip.py` from [here](https://bootstrap.pypa.io/get-pip.py)
3. Using the terminal, navigate to the directory with the `get-pip.py` file and run `python get-pip.py`
    - For macs, if you get a Permission Error, try `sudo python get-pip.py`
    - For windows users, try opening your command prompt/powershell window as Administrator, and reissue the `python get-pip.py` command.

## To install necessary files to run the program
- Prerequisites: Python 2.7, pip
1. Navigate to the directory with the `c_builder` directory in it.
2. Run `pip install -r requirements.txt` from the terminal.
    - If you get a Permission Error, try `pip install -r requirements.txt --user`

## To run the program
1. From the terminal, run the `run_server` script with `./run_server`
    - If on Windows use `run_server.ps1`, which can be run with `run_server.ps1`
2. Navigate to `localhost:5000` in your web browser.
