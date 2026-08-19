#!/bin/bash
##############################################################################
# Setup script for CEPCSW:
# - setup the external libraries
#
# Usage:
# $ source setup.sh
#
# Author: Tao Lin <lintao@ihep.ac.cn>
##############################################################################

THISSCRITDIR=$(dirname $(readlink -e "${BASH_SOURCE[0]}" 2>/dev/null) 2>/dev/null) # Darwin readlink doesnt accept -e

function info:() {
    echo "INFO: $*" 1>&2
}

function error:() {
    echo "ERROR: $*" 1>&2
}

function lcg-version() {
    echo ${CEPCSW_LCG_VERSION}
}

function cepcsw-base() {
    echo /cvmfs/cepcsw.ihep.ac.cn/prototype
}

function cepcsw-external() {
    echo $(cepcsw-base)/releases/externals/$(lcg-version)
}

function setup-external() {
    local f

    # try new os
    f=$(cepcsw-external)/setup-${CEPCSW_LCG_VERSION}-${CEPCSW_LCG_OS}.sh
    if [ -f "$f" ]; then
        info: "Setup CEPCSW externals (os: $CEPCSW_LCG_OS, lcg: $CEPCSW_LCG_VERSION): $f"
        source $f
        return
    fi

    # legacy
    if [ "$CEPCSW_LCG_OS" = "centos7" ]; then
        f=$(cepcsw-external)/setup-${CEPCSW_LCG_VERSION}.sh
        if [ -f "$f" ]; then
            info: "Setup CEPCSW externals (os: $CEPCSW_LCG_OS, lcg: $CEPCSW_LCG_VERSION): $f"
            source $f
            return
        fi
    fi

    error: "Failed to find setup script for CEPCSW_LCG_VERSION: ${CEPCSW_LCG_VERSION}, CEPCSW_LCG_OS: ${CEPCSW_LCG_OS}"

}

function setup-install-area() { 
    local installarea=$THISSCRITDIR/InstallArea/$CEPCSW_LCG_PLATFORM
    if [ ! -d "$installarea" ]; then
        info: "The InstallArea is not created yet. Please run ./build.sh first, then source setup.sh again."
        return 0
    fi

    export PATH=$installarea/bin:$PATH
    export LD_LIBRARY_PATH=$installarea/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$installarea/lib:$PYTHONPATH
    export PYTHONPATH=$installarea/python:$PYTHONPATH
    export ROOT_INCLUDE_PATH=$installarea/include:$ROOT_INCLUDE_PATH

    local extrasetupscript=$installarea/setup.sh
    if [ -f "$extrasetupscript" ]; then
        source $extrasetupscript
    fi

    info: "Setup CEPCSW: $installarea"
}

function detect-os() {
    local os=$(source /etc/os-release; echo $ID)
    local os_ver=$(source /etc/os-release; echo $VERSION_ID)

    local cepcsw_lcg_os
    
    case "$os" in
        almalinux)
            cepcsw_lcg_os=el
        ;;
        centos)
            cepcsw_lcg_os=centos
        ;;
    esac

    case "$os_ver" in
        7*)
            cepcsw_lcg_os=${cepcsw_lcg_os}7
            ;;
        9*)
            cepcsw_lcg_os=${cepcsw_lcg_os}9
            ;;
    esac

    export CEPCSW_LCG_OS=${cepcsw_lcg_os}

}

##############################################################################
# Parse the command line options
##############################################################################


# CEPCSW_LCG_OS: centos7 or el9
#if [ -z "$CEPCSW_LCG_OS" ]; then
#    # try to detect os
#    detect-os
#fi
#
## if user don't specify the version, try to detect according to the current os
#if [ -z "$CEPCSW_LCG_VERSION" ]; then
#    case "$CEPCSW_LCG_OS" in
#        el9)
#            CEPCSW_LCG_VERSION=105.0.0
#            ;;
#        *)
#            CEPCSW_LCG_VERSION=103.0.2
#            ;;
#    esac
#fi
CEPCSW_LCG_OS=el9
CEPCSW_LCG_VERSION=105.0.0
export CEPCSW_LCG_VERSION


setup-external
setup-install-area
