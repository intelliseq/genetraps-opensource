#!/bin/bash

# how to use
# ./save-secret.sh KEY VALUE bucket
KEY=$1
VALUE=${!KEY}
BUCKET=param-store
echo $KEY
echo $VALUE
echo $BUCKET
aws kms encrypt --key-id alias/param-store-key --plaintext $VALUE --query CiphertextBlob --output text --profile genetraps | base64 -d | aws s3 cp - s3://$BUCKET/$KEY --profile genetraps

#aws kms decrypt --ciphertext-blob fileb://<(aws s3 cp s3://kms-store/$KEY -) --output text --query Plaintext | base64 -d
