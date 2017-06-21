#!/bin/bash

docker rmi -f jorvis/biocode
docker build --no-cache -t jorvis/biocode .
docker images

echo "Now run: docker tag <newest tag here> jorvis/biocode:latest"
echo "         docker push jorvis/biocode"
