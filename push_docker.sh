#!/bin/sh

set -ev

VERSION=`cat VERSION.txt`

docker push yanay/combined_benchmark:${VERSION}
docker push yanay/combined_benchmark:latest
