
## Setup
Install [docker](https://www.docker.com/). On Windows/OSX, you may need to modify preferences of the Docker service to utilize more than 2 CPU cores and allow up to 8GB RAM usage. As necessary, modify `numCores <- 6` in `loadData.R` to specify the number of cores to run during the forecasting.

Build the docker image
```
docker build . -t m4:5.0 --build-arg BUILD_DATE=2018-04-11 --compress
```

Note that the `BUILD_DATE` should _not_ be modified if a reproducible build is desired since this pins the versions of the packages installed to those available at this date.

## Forecasts

Creating the forecasts is computationally expensive and could take around one week on a fast multi-core system. Run the image:
```
docker run -it --name forecastHybrid m4:5.0
```
If the process completes successfully, you will receive a message in the console:
```
Point forecasts and prediction intervals created successfully!
Pressing enter will exit and erase all results.
```
Before exiting, you will wish to copy the results from the Docker container to your filesystem. The forecasts will appear **inside** the docker container in the `mean/` directory; the upper/lower 95% prediction intervals will appear in `upper/` and `lower/`, respectively. In another terminal (_not_ inside the docker container), run the script `./copyFromContainer.sh` to copy the forecasts from the container to your filesystem. The `mean/`, `upper/', and `lower/` directories of your local filesystem will now contain the forecasts and prediction intervals, and you can safely exit and destroy the container.
