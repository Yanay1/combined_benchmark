#!/bin/sh

set -ev

VERSION=`cat VERSION.txt`

docker build -t yanay/combined_benchmark:$VERSION .
docker build -t yanay/combined_benchmark:latest .

