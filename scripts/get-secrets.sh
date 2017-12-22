#!/bin/bash
#
# use . ./get-secrets.sh kms-store 

BUCKET=$1
echo "Retrieving keys"

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     o=d;;
    Darwin*)    o=D;;
esac


for KEY in `aws s3 ls s3://$BUCKET | tr -s ' ' | cut -d " " -f 4`
do
    echo $KEY
    VALUE=`aws kms decrypt --ciphertext-blob fileb://<(aws s3 cp s3://kms-store/$KEY -) --output text --query Plaintext | base64 -$o`
    eval "export $KEY=$VALUE"
done
