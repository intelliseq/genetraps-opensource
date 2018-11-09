#!/bin/bash
#
# use . source encrypt.sh directory/file

PASSED=$1
DESTINATION="secret-data"
PASS=$SECRET_PASS
NAME=`echo "$PASSED" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`

if [[ -d $PASSED ]]; then
  WHERE=`pwd`
  cd $PASSED/..
  zip -r $NAME $NAME
else
  zip -j $NAME $PASSED
fi

ENCARCHIVE="$NAME.zip.enc"
echo $PASS | openssl enc -in "$NAME.zip" -aes-256-cbc -pass stdin > $ENCARCHIVE

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. >/dev/null && pwd )"
if [ ! -d "$DIR/$DESTINATION" ]; then
  mkdir $DIR/$DESTINATION
fi
mv $ENCARCHIVE $DIR/$DESTINATION
rm "$NAME.zip"

if [[ -d $PASSED ]]; then
  cd $WHERE
fi
