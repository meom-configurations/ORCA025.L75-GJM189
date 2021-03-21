#!/bin/ksh
chkdir () { if [ ! -d $1 ] ; then mkdir $1 ; fi ; }

  mkdir -p /gpfsstore/rech/cli/rcli002
  chkdir ORCA025.L75
  chkdir ORCA025.L75/ORCA025.L75-I
  chkdir ORCA025.L75/ORCA025.L75-GJM189-S
  chkdir ORCA025.L75/ORCA025.L75-GJM189-R
  chkdir ORCA025.L75/ORCA025.L75-GJM189-MEAN
