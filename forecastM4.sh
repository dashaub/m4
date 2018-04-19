#!/bin/bash

set -euo pipefail

Rscript forecastM4.R
echo 'Point forecasts and prediction intervals created successfully!'
echo 'Pressing enter will exit the container and erase all results.'
read
exit
