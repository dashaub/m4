Build the docker image

```
docker build . -t m4:4.0 --build-arg BUILD_DATE=2018-04-10 --compress
```

Run the image
```
docker run -it --name forecastHybrid m4:4.0
```
