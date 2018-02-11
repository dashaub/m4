#!/bin/bash
set -e
mkdir -p Data
cd Data/
wget https://www.m4.unic.ac.cy/wp-content/uploads/2017/12/M4DataSet.rar
unrar x M4DataSet.rar
rm M4DataSet.rar
