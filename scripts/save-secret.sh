#!/bin/bash

# how to use
# source save-secret.sh KEY
KEY=$1
VALUE=${!KEY}
BUCKET=param-store
echo $KEY
#echo $VALUE
echo $BUCKET

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     o=d;;
    Darwin*)    o=D;;
esac

aws kms encrypt --key-id alias/param-store-key --plaintext $VALUE --query CiphertextBlob --output text --profile genetraps | base64 -$o | aws s3 cp - s3://$BUCKET/$KEY --profile genetraps

#aws kms decrypt --ciphertext-blob fileb://<(aws s3 cp s3://kms-store/$KEY -) --output text --query Plaintext | base64 -d
