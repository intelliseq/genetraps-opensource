#!/bin/bash

# how to use
# ./save-file.sh KEY(name-in-bucket) VALUE(file)
KEY=$1
#VALUE=$2
BUCKET=fil-store
echo $KEY
#echo $VALUE
echo $BUCKET


aws s3 cp $2 s3://$BUCKET/$KEY --profile genetraps
