Build the docker image

```
docker build . -t m4:3.0 --build-arg BUILD_DATE=2018-04-09 --compress
```

Run the image
```
docker run -it --name forecastHybrid m4:3.0
```
