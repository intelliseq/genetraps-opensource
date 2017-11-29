echo "Running genetraps"
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure --profile kms
. scripts/get-secrets.sh
dx login --token $DNANEXUS_TOKEN --noprojects
echo `dx ls`
exit 0
