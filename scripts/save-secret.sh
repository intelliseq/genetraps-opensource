#!/bin/bash

# how to use
# ./save-secret.sh KEY VALUE bucket
KEY=$1
VALUE=$2
BUCKET=$3
echo $KEY
echo $VALUE
echo $BUCKET
aws kms encrypt --key-id alias/stash --plaintext $VALUE --query CiphertextBlob --output text | base64 -d | aws s3 cp - s3://$BUCKET/$KEY

#aws kms decrypt --ciphertext-blob fileb://<(aws s3 cp s3://kms-store/$KEY -) --output text --query Plaintext | base64 -d
