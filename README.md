# 8-ID-I XPCS Bluesky Instrument

## Installation Steps
*aps_8id_bs_instrument* can also use *conda* for dependency management, and
*setuptools* for installation and development.

First, download the package from github:

```bash
git clone https://github.com/aps-8id-dys/bluesky
cd bluesky
```

```bash
export ENV_NAME=bs_8id_main
conda create -y -n $ENV_NAME python=3.11 pyepics apsu::aps-dm-api
conda activate $ENV_NAME
pip install -e ."[all]"
```

## IPython console
To start the bluesky instrument session in a ipython execute the below command in a terminal
```bash
ipython
```
Inside the ipython console execute
```py
from instrument.startup import *
```
## Jupyter notebook

Start JupyterLab, a Jupyter notebook server, or a notebook, VSCode.

Start the data acquisition:

```py
from instrument.startup import *
```

## Sim Plan Demo
To run some simulated plans that ensure the installation worked as expected please run the below inside an ipython session or a jupyter notebook after starting the data acquisition
```py
RE(sim_print_plan())
RE(sim_count_plan())
RE(sim_rel_scan_plan())
```

See this [example](./docs/source/demo.ipynb).

## Configuration files
The files that can be configured to adhere to your preferences are:
- `configs/iconfig.yml` - configuration for data collection
- `configs/logging.yml` - configuration for session logging to console and/or files
- `qs/qs-config.yml`    - contains all configuration of the QS host process. See the [documentation](https://blueskyproject.io/bluesky-queueserver/manager_config.html) for more details of the configuration.

## queueserver

The queueserver has a host process that manages a RunEngine. Client sessions
will interact with that host process.

### Run a queueserver host process

Use the queueserver host management script to start the QS host process.  The below option stops the server (if it
is running) and then starts it.  This is the usual way to (re)start the QS host
process. Using the below command the process runs in the background.

```bash
./qs/qs_host.sh restart
```

### Run a queueserver client GUI
To run the gui client for the queueserver you can use the below command inside the terminal
```bash
queue-monitor &
```

### Shell script explained

A [shell script](./qs/qs_host.sh) is used to start the QS host process. Below are all the command options, and what they do.
```bash
(bstest) $ ./qs/qs_host.sh help
Usage: qs_host.sh {start|stop|restart|status|checkup|console|run} [NAME]

    COMMANDS
        console   attach to process console if process is running in screen
        checkup   check that process is running, restart if not
        restart   restart process
        run       run process in console (not screen)
        start     start process
        status    report if process is running
        stop      stop process

    OPTIONAL TERMS
        NAME      name of process (default: bluesky_queueserver-)
```

Alternatively, run the QS host's startup comand directly within the `qs/`
subdirectory.

```bash
cd ./qs
start-re-manager --config=./qs-config.yml
```

## Testing

Use this command to run the test suite locally:
```bash
pytest -vvv --lf ./src
```

# Warnings
##  For the Bluesky Queueserver.

The QS host process writes files into the qs directory. This directory can be
relocated. However, it should not be moved into the instrument package since
that might be installed into a read-only directory.


## OLD README (review)

### Sanity check test script
The below script will allow you to check if you are on the private subnet as
well as check if APS Data Management tools are installed properly.
```bash
python3 ./scripts/user/check_environment_test.py
```
## Running Bluesky Session
### With Ipython

```bash
./scripts/bs_ipy_starter.sh
```

Then Inside the ipython shell

```bash
RE(demo_sim_1d())
```

### With Queserver

Inside one terminal

```bash
./scripts/bs_qs_screen_starter.sh run
```

Inside another terminal

```bash
qserver environment open
```

```bash
qserver queue add plan '{"name": "demo_sim_1d"}'
qserver queue start
```

## Useful Bluesky Commands

command | description
--- | ---
`listObjects()` | show all ophyd devices available
`device_name.component_names` | show all parts of `device_name` (such as motors)
`device_name.summary()` | details of `device_name`

## SPEC to Bluesky Cheatsheet
Please refer to the cheatsheet below in case you need or desire to run commands directly through the ipython session
https://bcda-aps.github.io/bluesky_training/howto/bluesky_cheat_sheet.html
