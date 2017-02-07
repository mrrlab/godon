#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
buildstamp=$(date -u '+%Y-%m-%d_%H:%M:%S')
githash=$(git -C $DIR rev-parse HEAD)
gitbranch=$(git -C $DIR rev-parse --abbrev-ref HEAD)
go install -ldflags "--extldflags=-static -X main.buildstamp=$buildstamp -X main.githash=$githash -X main.gitbranch=$gitbranch" bitbucket.org/Davydov/godon/godon
