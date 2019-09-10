#!/bin/bash
#
# use . source decrypt.sh name.enc

PASSED=$1
DESTINATION="secret-data/decrypted"
PASS=$SECRET_PASS

DECARCHIVE=`echo "${PASSED%".enc"}" | rev | cut -d'/' -f1 | rev`
echo $PASS | openssl enc -in $PASSED -d -aes-256-cbc -pass stdin > $DECARCHIVE

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. >/dev/null && pwd )"
if [ ! -d "$DIR/$DESTINATION" ]; then
  mkdir -p $DIR/$DESTINATION
fi

yes | unzip $DECARCHIVE -d $DIR/$DESTINATION
rm $DECARCHIVE
