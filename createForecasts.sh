#!/bin/bash

set -euo pipefail

Rscript loadData.R
echo 'Point forecasts and prediction intervals created successfully!'
echo 'Pressing enter will exit and erase all results.'
read
exit
