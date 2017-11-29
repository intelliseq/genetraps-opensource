#!/bin/bash
#
# use . ./get-secrets.sh kms-store 

BUCKET=$1

for KEY in `aws s3 ls s3://$BUCKET | tr -s ' ' | cut -d " " -f 4`
do
    VALUE=`aws kms decrypt --ciphertext-blob fileb://<(aws s3 cp s3://kms-store/$KEY -) --output text --query Plaintext | base64 -d`
    eval "$KEY=$VALUE"
done
