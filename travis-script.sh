LOG_PREFIX="===== TRAVIS LOG ===== "
LOG_APP=" no-app: "

echo $LOG_PREFIX"Starting travis script"
source scripts/get-secrets.sh kms-store
echo $LOG_PREFIX"Checking first five letters of token: "`echo $DNANEXUS_TOKEN | cut -c1-5`
printf "0\n" | dx login --token $DNANEXUS_TOKEN
#echo `dx ls`

### checking for errors
function check {
if [ $? -eq 0 ]; then
    echo $LOG_PREFIX"service up"
else
    echo $LOG_PREFIX"service failed"
    exit 1
fi
} 

### BUILDING API_DX ###
LOG_APP="api-dx: "
echo $LOG_PREFIX $LOG_APP "building..."
gradle build docker -p api-dx/
echo $LOG_PREFIX $LOG_APP "running docker..."
docker run -d -p 8080:8080 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN -t pl.intelliseq.genetraps.api.dx/api-dx:latest
echo $LOG_PREFIX $LOG_APP "waiting for service..."
./scripts/wait-for-service.sh localhost:8080/hello 60
echo $LOG_PREFIX $LOG_APP "checking for error..."
check
echo $LOG_PREFIX $LOG_APP "testing..."
echo $(curl localhost:8080/touch)

### BUILDING CLINET_EXPLORARE ###
LOG_APP="client-explorare: "
echo $LOG_PREFIX $LOG_APP "compiling..."
mkdir ./client-explorare/dist
npm --prefix ./client-explorare install ./client-explorare
echo $LOG_PREFIX $LOG_APP "testing..."
npm --prefix client-explorare run test
echo $LOG_PREFIX $LOG_APP "client-explorare: dockerizing..."
docker build -t pl.intelliseq.genetraps.client.explorare/client-explorare:latest client-explorare
echo $LOG_PREFIX $LOG_APP "client-explorare: running docker..."
docker run -d -p 8081:8081 pl.intelliseq.genetraps.client.explorare/client-explorare:latest
echo $LOG_PREFIX $LOG_APP "client-explorare: waiting for service..."
./scripts/wait-for-service.sh localhost:8081/hello 60
echo $LOG_PREFIX $LOG_APP "client-explorare: checking for error..."
check


echo $LOG_PREFIX "pushing to AWS"
ecs-cli configure profile --profile-name genetraps --access-key $ECS_CLI_GENETRAPS_DEV_KEY_ID --secret-key $ECS_CLI_GENETRAPS_DEV_ACCESS_KEY
ecs-cli configure --cluster GENETRAPS-DEV --default-launch-type FARGATE --region us-east-1 --config-name ECS_CLI_GENETRAPS_DEV_CONF
mkdir ~/.aws/
printf "[genetraps]\naws_access_key_id="$ECS_CLI_GENETRAPS_DEV_KEY_ID"\naws_secret_access_key="$ECS_CLI_GENETRAPS_DEV_ACCESS_KEY >> ~/.aws/credentials
printf "[profile genetraps]\nregion="$AWS_REGION"\noutput=json" >> ~/.aws/config
echo "ecs cli"
echo `cat ~/.aws/credentials | cut -c -26`
echo `cat ~/.aws/credentials | wc -l`
echo `ls -la ~/.aws/`
aws ecr get-login --profile genetraps | echo
`aws ecr get-login --profile genetraps`
echo `ecs-cli ps --ecs-profile genetraps`
echo `ecs-cli ps --aws-profile genetraps`

TAG_API_DX=`ls -alR --full-time api-dx/ -Ibin -Ibuild -I.* | sha1sum | cut -d " " -f1`
TAG_CLIENT_EXPLORARE=`ls -alR --full-time client-explorare/ -Inode_modules -Idist -Ietc -I.* | sha1sum | cut -d " " -f1`
#echo $TAG_API_DX
#echo $TAG_CLIENT_EXPLORARE
IMAGE_API_DX=pl.intelliseq.genetraps.api.dx/api-dx:latest
IMAGE_CLIENT_EXPLORARE=pl.intelliseq.genetraps.client.explorare/client-explorare:latest
#echo $IMAGE_CLIENT_EXPLORARE
TAG_API_DX=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-api-dx:"$TAG_API_DX
TAG_CLIENT_EXPLORARE=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-client-explorare:"$TAG_CLIENT_EXPLORARE
docker tag $IMAGE_API_DX $TAG_API_DX
docker tag $IMAGE_CLIENT_EXPLORARE $TAG_CLIENT_EXPLORARE
ecs-cli push $TAG_API_DX --ecs-profile genetraps
ecs-cli push $TAG_CLIENT_EXPLORARE --ecs-profile genetraps
#echo $TAG_CLIENT_EXPLORARE 
cat aws-conf/docker-compose-template.yml | sed 's@clientExplorareImageTag@'"$TAG_CLIENT_EXPLORARE"'@' > docker-compose.yml
ecs-cli compose --project-name genetraps-client-explorare -f docker-compose.yml --ecs-params ./aws-conf/ecs-params.yml service up --target-group-arn "arn:aws:elasticloadbalancing:us-east-1:"$AWS_ACCOUNT_ID":targetgroup/genetraps-explorare-client/"$AWS_GENETRAPS_TARGET_GROUP --container-name client-explorare --container-port 8081 --ecs-profile genetraps

#./scripts/update-repo.sh

exit 0
