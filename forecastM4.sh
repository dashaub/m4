#!/bin/bash

set -euo pipefail

cd ~/m4
~/m4/DownloadData.sh

Rscript /root/m4/forecastM4.R $1
echo 'Point forecasts and prediction intervals created successfully!'
echo 'Pressing enter will exit the container and erase all results.'
read
exit
