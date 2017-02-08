#!/bin/bash
function print_help {
	echo "Usage: $0 [-d] [-f]" >&2
	echo >&2
	echo "-d compiles a dynamic binary" >&2
	echo "-f force recompilation by deleting the binary first" >&2
}

function godepinstalled {
	if ! go list "$1" &> /dev/null || ! test -f $(go list -f '{{.Target}}'  "$1")
	then
		echo You need to install $1 first >&2
		exit
	fi
}

function goisold {
	ver=$(go version | awk '{print $3}')
	echo $ver | grep -q -E 'go1\.[0-6](\.[0-9]*)?$'
}

STATIC="--extldflags=-static"

while getopts "dhf" opt
do
  case $opt in
    d)
      echo "building a dynamic binary" >&2
      STATIC=""
      ;;
    f)
      FORCE=1
      ;;
    h)
      print_help
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      print_help
      exit
      ;;
  esac
done

godepinstalled github.com/gonum/blas/cgo

# prior to go1.7 go-lbfgsb requires a manual installation
goisold && \
	godepinstalled github.com/idavydov/go-lbfgsb

PKG_NAME=bitbucket.org/Davydov/godon/godon
BINARY_NAME=$(basename $PKG_NAME)
BINARY_PATH=$GOPATH/bin/$BINARY_NAME
go get -d $PKG_NAME
GOPATH=${GOPATH-$HOME/go}
DIR="$GOPATH/src/$PKG_NAME"
buildstamp=$(date -u '+%Y-%m-%d_%H:%M:%S')
githash=$(git -C "$DIR" rev-parse HEAD)
gitbranch=$(git -C "$DIR" rev-parse --abbrev-ref HEAD)
if [ -n "$FORCE" ]
then
	rm -f $BINARY_PATH
fi
go get -ldflags "$STATIC -X main.buildstamp=$buildstamp -X main.githash=$githash -X main.gitbranch=$gitbranch" $PKG_NAME && \
	echo Installed $BINARY_NAME to $BINARY_PATH >&2
