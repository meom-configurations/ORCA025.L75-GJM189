#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
#
# ==========
# agrifpp.sh
# ==========
#
# ----------------------------
# Preform AGrif pre-processing
# ----------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ agrifpp.sh
#
#
# DESCRIPTION
# ===========
#
#
# Preprocess file using the conv in OPAFILES directory
# Standard preprocessed files are stored in OPAFILES/ppsrc/nemo
# Source files are stored under OPAFILES/obj
# Include filess  in OPAFILES/inc
# Note that agrif2model.F90 should not be preprocess (standard one) 
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./agrifpp.sh FILE_TO_PROCESS
# 
# TODO
# ====
#
# option debug
#
#
# EVOLUTIONS
# ==========
#
# $Id: agrifpp.sh 2143 2010-10-04 12:49:55Z rblod $
#
#
#
#   * creation
#
#-
MYFILE=$(basename "$1")
if [ "$MYFILE" == "agrif2model.f90" ];then
   \cp ${NEMO_TDIR}/${NEW_CONF}/WORK/${MYFILE/.f90/.F90} ${NEMO_TDIR}/${NEW_CONF}/OPAFILES/obj/$MYFILE
else
cd ${NEMO_TDIR}/${NEW_CONF}/OPAFILES/ppsrc/nemo ; ${NEMO_TDIR}/${NEW_CONF}/OPAFILES/conv ${NEMO_TDIR}/${NEW_CONF}/OPAFILES/agrif_opa.in -rm -incdir ${NEMO_TDIR}/${NEW_CONF}/OPAFILES/inc -comdirout ${NEMO_TDIR}/${NEW_CONF}/OPAFILES/obj -convfile ${MYFILE} > /dev/null 
fi