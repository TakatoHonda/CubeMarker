#!/bin/sh

#
# demo.sh --- written by Takato Honda
# 

make clean
make
INPUT="./_list/sample.txt"
OUTDIR="./_out/"

#----------------------#
echo "======================================"
echo " _      _   _            _      _  _  "
echo "/  | | |_) |_ |\/|  /\  |_) |/ |_ |_) "
echo "\_ |_| |_) |_ |  | /--\ | \ |\ |_ | \ "
echo "                                      "
echo "            --------------------      "
echo "         / |  /\  _  /\_  /\  _/|     "
echo "        -------------------- /  |     "
echo "     / |  /\  _  /\_  /\  _/|---      "
echo "    -------------------- /  | /       "
echo "   |  /\  _  /\_  /\  _/|---          "
echo "   | /  \/ \/   \/  \/  | /           "
echo "    --------------------              "
echo "======================================"
#----------------------#

./cubemarker $INPUT $OUTDIR