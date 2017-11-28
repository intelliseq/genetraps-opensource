echo "Running genetraps"
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure --profile kms
echo `aws kms list-aliases`
exit 0
