#!/bin/bash
#
# use . ./get-jwt.sh

KEY=JWT
BUCKET=fil-store
echo "Retrieving jwt.jks file and creating public.cert"
echo $KEY

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. >/dev/null && pwd )"
cd $DIR
#echo $DIR

aws s3 cp s3://$BUCKET/$KEY $DIR/api-dx/src/main/resources/jwt.jks --profile genetraps
cp $DIR/api-dx/src/main/resources/jwt.jks $DIR/api-security/src/main/resources/jwt.jks

echo $SECRET_KEY | keytool -list -rfc --keystore $DIR/api-dx/src/main/resources/jwt.jks | openssl x509 -inform pem -pubkey | head -n 9 | cat > $DIR/api-dx/src/main/resources/public.cert
cp $DIR/api-dx/src/main/resources/public.cert $DIR/api-security/src/main/resources/public.cert
