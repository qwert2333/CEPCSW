#!/bin/bash
##############################################################################
# Run script for CEPCSW:
# - run a simple job
#
# Usage:
# $ ./run.sh Examples/options/helloalg.py
# or:
# $ 
#
# Author: Tao Lin <lintao@ihep.ac.cn>
##############################################################################

function info:() {
    echo "INFO: $*" 1>&2
}

function error:() {
    echo "ERROR: $*" 1>&2
}

function check-cepcsw-envvar() {
    # CEPCSWEXTERNAL is defined in /cvmfs/cepcsw.ihep.ac.cn/prototype/releases/externals/
    if [ -z "${CEPCSWEXTERNAL}" ]; then
        error: "The CEPCSW is not setup. Please source setup.sh."
        return 1
    fi
}

function build-dir() {
    local blddir=$WorkDIR/build

    echo $blddir
}

function check-working-builddir() {
    local blddir=$(build-dir)
    if [ ! -d "$blddir" ]; then
        mkdir $blddir || {
            error: "Failed to create $blddir"
            return 1
        }
    fi
}

function run-job() {
    local blddir=$(build-dir)

    $blddir/run gaudirun.py $*
}

##############################################################################
# Parse the command line options
##############################################################################

# The current default platform
lcg_platform=x86_64-centos7-gcc8-opt
<<<<<<< HEAD:run/run.sh
lcg_version=98.0.0
=======
lcg_version=101.0.1

bldtool=${CEPCSW_BLDTOOL} # make, ninja # set in env var
>>>>>>> d084d9f06570119e90ad0e4714585dc2f4f69f1a:run.sh

check-cepcsw-envvar || exit -1

check-working-builddir || exit -1

run-job $* || exit -1
