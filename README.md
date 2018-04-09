Build the docker image

```
docker build . -t m4:1.0 --build-arg BUILD_DATE=2018-04-08
```

Run the image
```
docker run -it m4:1.0
```
