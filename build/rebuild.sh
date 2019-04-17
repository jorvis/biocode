#!/bin/bash

# This should match the PyPi version, which needs to be available first.
IMAGE_VERSION='0.9.0'

docker rmi -f jorvis/biocode
docker build --no-cache -t jorvis/biocode:latest -t jorvis/biocode:${IMAGE_VERSION} .
docker images

echo "If ready for release, run: "
echo "  docker push jorvis/biocode:latest"
echo "  docker push jorvis/biocode:${IMAGE_VERSION}"
