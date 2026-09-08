#!/bin/bash
#
# Alicat PCD pressure-controller IOC (8idAlicat:PCD1 and 8idAlicat:PCD2).
#
# Wrapper around the IOC's own 8idAlicat.pl so you do not have to remember, or
# cd into, the path below. Every argument is passed straight through, so any
# command 8idAlicat.pl grows in future works here too without editing this file.
#
#   start_Alicat.sh status     is it running?
#   start_Alicat.sh start      start it
#   start_Alicat.sh stop       stop it
#   start_Alicat.sh restart    stop then start
#   start_Alicat.sh caqtdm     open the two Alicat_PCD.ui control screens
#   start_Alicat.sh console    attach to the IOC shell (screen; ctrl-a d to detach)
#   start_Alicat.sh            with no argument, prints the IOC's own usage
#
# The IOC talks to the controllers over a MOXA terminal server at
# 10.54.116.100 ports 4001 (PCD1) and 4002 (PCD2), so it must run somewhere
# that can reach that address.

IOC_DIR="${ALICAT_IOC_DIR:-/net/s8iddserv/xorApps/epics/synApps_6_3/ioc/8idAlicat/iocBoot/ioc8idAlicat/softioc}"
IOC_SCRIPT="8idAlicat.pl"

if [ ! -d "$IOC_DIR" ]; then
    echo "start_Alicat.sh: IOC directory not found:" >&2
    echo "    $IOC_DIR" >&2
    echo "Is /net/s8iddserv mounted here? Override with ALICAT_IOC_DIR=<path>." >&2
    exit 1
fi

if [ ! -x "$IOC_DIR/$IOC_SCRIPT" ]; then
    echo "start_Alicat.sh: $IOC_DIR/$IOC_SCRIPT is missing or not executable." >&2
    exit 1
fi

# cd first: 8idAlicat.pl resolves its own commands/ directory relative to cwd.
cd "$IOC_DIR" || exit 1
exec ./"$IOC_SCRIPT" "$@"
