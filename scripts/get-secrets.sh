#!/bin/bash
#
# use . ./get-secrets.sh kms-store 

BUCKET=param-store
echo "Retrieving keys"

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     o=d;;
    Darwin*)    o=D;;
esac


for KEY in `aws s3 ls s3://$BUCKET --profile genetraps | tr -s ' ' | cut -d " " -f 4`
do
    echo $KEY
    VALUE=`aws kms decrypt --ciphertext-blob fileb://<(aws s3 cp s3://param-store/$KEY - --profile genetraps) --output text --query Plaintext --profile genetraps | base64 -$o`
    eval "export $KEY=$VALUE"
done
