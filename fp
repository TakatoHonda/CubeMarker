#!/bin/sh

[ -n "$1" ] && cd $1
#-v : for linux 
for F in `ls -v`
do
  echo `pwd`/$F
done

