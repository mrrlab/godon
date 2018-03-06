#!/bin/bash
GO=go

function print_help {
	echo "Usage: $0 [-d] [-f]" >&2
	echo >&2
	echo "-d compiles a dynamic binary" >&2
	echo "-f force recompilation by deleting the binary first" >&2
	exit
}

function godepinstalled {
	if ! $GO list "$1" &> /dev/null || ! test -f $($GO list -f '{{.Target}}'  "$1")
	then
		echo Package $1 is not installed >&2
		exit
	fi
}

function goisold {
	ver=$($GO version | awk '{print $3}')
	echo $ver | grep -q -E 'go1\.[0-6](\.[0-9]*)?$'
}

if [[ "$OSTYPE" == "linux-gnu" ]]
then
    STATIC="--extldflags=-static"
else
    ## Mac OS X doesn't support static binaries
    STATIC=""
fi

if [[ "$BLAS" == "" ]]
then
   if [[ "$OSTYPE" == "linux-gnu" ]]
   then
       # OpenBLAS by default for GNU/Linux
       BLAS="-lopenblas"
   elif [[ "$OSTYPE" == "darwin"* ]]
   then
       BLAS="-framework Accelerate"
   else
       echo Unknown platform, using openblas >&2
       BLAS="-lopenblas"
   fi
fi
   

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

CGO_LDFLAGS="$CGO_LDFLAGS $BLAS" $GO get github.com/gonum/blas/cgo
godepinstalled github.com/gonum/blas/cgo

# prior to go1.7 go-lbfgsb requires a manual installation
goisold && \
	godepinstalled github.com/idavydov/go-lbfgsb

GOPATH=${GOPATH:-$HOME/go}
PKG_NAME=bitbucket.org/Davydov/godon/godon
BINARY_NAME=$(basename $PKG_NAME)
BINARY_PATH=$GOPATH/bin/$BINARY_NAME
$GO get -d $PKG_NAME
DIR="$GOPATH/src/$PKG_NAME"
buildstamp=$(date -u '+%Y-%m-%d_%H:%M:%S')
githash=$(git -C "$DIR" rev-parse HEAD)
gitbranch=$(git -C "$DIR" rev-parse --abbrev-ref HEAD)
if [ -n "$FORCE" ]
then
	rm -f $BINARY_PATH
fi
$GO get -ldflags "$STATIC -X main.buildstamp=$buildstamp -X main.githash=${githash:-$BITBUCKET_COMMIT} -X main.gitbranch=${gitbranch:-$BITBUCKET_BRANCH}" $PKG_NAME && \
	echo Installed $BINARY_NAME to $BINARY_PATH >&2
