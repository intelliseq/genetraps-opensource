#!/bin/bash
set -x
TOOL_ID=$1
DXWDL_VERSION=0.74.1

# GET SCRIPT DIRECTORY
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# DXWDL #
if [ ! -f /tmp/dxWDL/dxWDL-${DXWDL_VERSION}.jar ]; then
    wget https://github.com/dnanexus/dxWDL/releases/download/${DXWDL_VERSION}/dxWDL-${DXWDL_VERSION}.jar --directory-prefix /tmp/dxWDL
fi
java -jar /tmp/dxWDL/dxWDL-${DXWDL_VERSION}.jar compile ${DIR}/${TOOL_ID}/${TOOL_ID}.wdl -f -extras ${DIR}/extras -project genetraps-test
