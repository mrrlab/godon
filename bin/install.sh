#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
go install -ldflags "-X main.buildstamp=`date -u '+%Y-%m-%d_%H:%M:%S'` -X main.githash=`git -C $DIR rev-parse HEAD` -X main.gitbranch=`git -C $DIR rev-parse --abbrev-ref HEAD`" bitbucket.org/Davydov/godon/godon
