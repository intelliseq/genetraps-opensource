echo "Running genetraps"
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure
#printf $AWS_ACCESS_KEY_ID"\n"$AWS_SECRET_ACCESS_KEY"\n"$AWS_REGION"\njson" | aws configure --profile kms
echo "Starting"
#echo `aws describe-my-user-profile`
echo "token"`echo $DNANEXUS_TOKEN | cut -c1-5`
source scripts/get-secrets.sh kms-store
echo "token"`echo $DNANEXUS_TOKEN | cut -c1-5`
dx login --token $DNANEXUS_TOKEN --noprojects
echo `dx ls`

### BUILDING API_DX ###
gradle build docker -p api-dx/
docker run -d -p 8080:8080 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN -t pl.intelliseq.genetraps.api.dx/api-dx:latest
./scripts/wait-for-service.sh localhost:8080/hello 60
if [ $? -eq 0 ]; then
    echo "api-dx service up"
else
    echo "api-dx service failed"
    exit 1
fi
echo $(curl localhost:8080/touch)

exit 0
