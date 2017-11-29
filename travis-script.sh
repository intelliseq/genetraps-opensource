echo "Running genetraps"
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure --profile kms
echo "Starting"
echo `aws describe-my-user-profile`
#echo "token"`echo $DNANEXUS_TOKEN | cut -c1-5`
#. scripts/get-secrets.sh kms-store
#echo "token"`echo $DNANEXUS_TOKEN | cut -c1-5`
#dx login --token $DNANEXUS_TOKEN --noprojects
echo `dx ls`
exit 0
