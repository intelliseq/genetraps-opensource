LOG_PREFIX="===== TRAVIS LOG ===== "
LOG_APP=" no-app: "

echo $LOG_PREFIX"Starting travis script"
source scripts/get-secrets.sh
echo $LOG_PREFIX"Checking first five letters of token: "`echo $DNANEXUS_TOKEN_TEST | cut -c1-5`
printf "0\n" | dx login --token $DNANEXUS_TOKEN_TEST


echo $LOG_PREFIX"Cleaning DNANEXUS test space"
dx rm -r -a /samples/
#ls dx-apps/ | xargs -I {} bash -c "dx build dx-apps/{}"
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

### LOGGING TO AWS DOCKER HUB
echo $LOG_PREFIX "logging to ecr..."
`aws ecr get-login --profile genetraps --no-include-email`

### BUILDING CLIENT_INDEX ###
LOG_APP="client-index: "
echo $LOG_PREFIX $LOG_APP "setting tag"
CLIENT_INDEX_CHECKSUM=`find client-index -type f -exec md5sum {} \; | sort -k 2 | md5sum | sed 's/  -//g'`
CLIENT_INDEX_TAG=$AWS_ACCOUNT_ID".dkr.ecr."$AWS_REGION".amazonaws.com/genetraps-client-index:"$CLIENT_INDEX_CHECKSUM
echo $LOG_PREFIX $LOG_APP $CLIENT_INDEX_TAG
docker build client-index/ -t $CLIENT_INDEX_TAG
docker push $CLIENT_INDEX_TAG
cat aws-conf/docker-compose-template.yml | sed 's@clientIndexImageTag@'"$CLIENT_INDEX_TAG"'@' > docker-compose.yml

### ECS ###
ecs-cli configure --cluster genetraps --default-launch-type FARGATE --region us-east-1 --config-name genetraps
ecs-cli compose --project-name genetraps-client-index -f docker-compose.yml --ecs-params ./aws-conf/ecs-params.yml service up --target-group-arn "arn:aws:elasticloadbalancing:"$AWS_REGION":"$AWS_ACCOUNTID":targetgroup/client-index-target-group/"$ECS_CLI_TG_CLIENT_INDEX --container-name client-index --container-port $ECS_CLI_PORT_CLIENT_INDEX --aws-profile genetraps

### BUILDING API_DX ###
#LOG_APP="api-dx: "
#echo $LOG_PREFIX $LOG_APP "building..."
#gradle build docker -p api-dx/
#echo $LOG_PREFIX $LOG_APP "running docker..."
#docker run -d -p 8086:8086 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN_TEST -t pl.intelliseq.genetraps.api.dx/api-dx:latest
#echo $LOG_PREFIX $LOG_APP "waiting for service..."
#./scripts/wait-for-service.sh localhost:8086/hello 60
#echo $LOG_PREFIX $LOG_APP "checking for error..."
#check
#echo $LOG_PREFIX $LOG_APP "testing..."
#echo $(curl localhost:8086/touch)

#exit 0

### BUILDING CLIENT_EXPLORARE ###
#LOG_APP="client-explorare: "
#echo $LOG_PREFIX $LOG_APP "compiling..."
#mkdir ./client-explorare/dist
#npm --prefix ./client-explorare install ./client-explorare
#echo $LOG_PREFIX $LOG_APP "testing..."
#npm --prefix client-explorare run test
#echo $LOG_PREFIX $LOG_APP "dockerizing..."
#docker build -t pl.intelliseq.genetraps.client.explorare/client-explorare:latest client-explorare
#echo $LOG_PREFIX $LOG_APP "running docker..."
#docker run -d -p 8081:8081 pl.intelliseq.genetraps.client.explorare/client-explorare:latest
#echo $LOG_PREFIX $LOG_APP "waiting for service..."
#./scripts/wait-for-service.sh localhost:8081/hello 60
#echo $LOG_PREFIX $LOG_APP "checking for error..."
#check

### BUILDING CLIENT_DNATOKEN ###
#LOG_APP="client-dnatoken: "
#echo $LOG_PREFIX $LOG_APP "compiling..."
#mkdir ./client-dnatoken/dist
#npm --prefix ./client-dnatoken install ./client-dnatoken
#echo $LOG_PREFIX $LOG_APP "testing..."
#npm --prefix client-dnatoken run test
#echo $LOG_PREFIX $LOG_APP "dockerizing..."
#docker build -t pl.intelliseq.genetraps.client.dnatoken/client-dnatoken:latest client-dnatoken
#echo $LOG_PREFIX $LOG_APP "running docker..."
#docker run -d -p 8083:8083 pl.intelliseq.genetraps.client.dnatoken/client-dnatoken:latest
#echo $LOG_PREFIX $LOG_APP "waiting for service..."
#./scripts/wait-for-service.sh localhost:8083/hello 60
#echo $LOG_PREFIX $LOG_APP "checking for error..."
#check

#echo $LOG_PREFIX "pushing to AWS"
#ecs-cli configure profile --profile-name genetraps --access-key $ECS_CLI_GENETRAPS_DEV_KEY_ID --secret-key $ECS_CLI_GENETRAPS_DEV_ACCESS_KEY
#ecs-cli configure --cluster GENETRAPS-DEV --default-launch-type FARGATE --region us-east-1 --config-name ECS_CLI_GENETRAPS_DEV_CONF
#mkdir ~/.aws/
#printf "[genetraps]\naws_access_key_id="$ECS_CLI_GENETRAPS_DEV_KEY_ID"\naws_secret_access_key="$ECS_CLI_GENETRAPS_DEV_ACCESS_KEY >> ~/.aws/credentials
#printf "[profile genetraps]\nregion="$AWS_REGION"\noutput=json" >> ~/.aws/config
#echo "ecs cli"
#echo `cat ~/.aws/credentials | cut -c -26`
#echo `cat ~/.aws/credentials | wc -l`
#echo `ls -la ~/.aws/`
#aws ecr get-login --profile genetraps | echo
#`aws ecr get-login --profile genetraps`
#echo `ecs-cli ps --ecs-profile genetraps`
#echo `ecs-cli ps --aws-profile genetraps`

#TAG_API_DX=`ls -alR --full-time api-dx/ -Ibin -Ibuild -I.* | sha1sum | cut -d " " -f1`
#TAG_CLIENT_EXPLORARE=`ls -alR --full-time client-explorare/ -Inode_modules -Idist -Ietc -I.* | tr -s ' ' | cut -f1-5,9 | sha1sum | cut -d " " -f1`
#TAG_CLIENT_DNATOKEN=`ls -alR --full-time client-dnatoken/ -Inode_modules -Idist -Ietc -I.* | tr -s ' ' | cut -f1-5,9 | sha1sum | cut -d " " -f1`
#echo $TAG_API_DX
#echo $TAG_CLIENT_EXPLORARE
#IMAGE_API_DX=pl.intelliseq.genetraps.api.dx/api-dx:latest
#IMAGE_CLIENT_EXPLORARE=pl.intelliseq.genetraps.client.explorare/client-explorare:latest
#IMAGE_CLIENT_DNATOKEN=pl.intelliseq.genetraps.client.dnatoken/client-dnatoken:latest
#echo $IMAGE_CLIENT_EXPLORARE
#TAG_API_DX=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-api-dx:"$TAG_API_DX
#TAG_CLIENT_EXPLORARE=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-client-explorare:"$TAG_CLIENT_EXPLORARE
#TAG_CLIENT_DNATOKEN=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-client-dnatoken:"$TAG_CLIENT_DNATOKEN
#docker tag $IMAGE_API_DX $TAG_API_DX
#docker tag $IMAGE_CLIENT_EXPLORARE $TAG_CLIENT_EXPLORARE
#docker tag $IMAGE_CLIENT_DNATOKEN $TAG_CLIENT_DNATOKEN
#ecs-cli push $TAG_API_DX --ecs-profile genetraps
#ecs-cli push $TAG_CLIENT_EXPLORARE --ecs-profile genetraps
#ecs-cli push $TAG_CLIENT_DNATOKEN --ecs-profile genetraps
#echo $TAG_CLIENT_EXPLORARE
#cat aws-conf/docker-compose-template.yml | sed 's@clientExplorareImageTag@'"$TAG_CLIENT_EXPLORARE"'@' | sed 's@clientDnatokenImageTag@'"$TAG_CLIENT_DNATOKEN"'@' > docker-compose.yml
#ecs-cli compose --project-name genetraps-client-explorare -f docker-compose.yml --ecs-params ./aws-conf/ecs-params.yml service up --target-group-arn "arn:aws:elasticloadbalancing:us-east-1:"$AWS_ACCOUNT_ID":targetgroup/genetraps-explorare-client/"$AWS_CLIENT_EXPLORAE_TARGET_GROUP --container-name client-explorare --container-port 8081 --ecs-profile genetraps
#cat aws-conf/docker-compose-template.yml | sed 's@clientDnatokenImageTag@'"$TAG_CLIENT_DNATOKEN"'@' > docker-compose.yml
#ecs-cli compose --project-name genetraps-client-dnatoken -f docker-compose.yml --ecs-params ./aws-conf/ecs-params.yml service up --target-group-arn "arn:aws:elasticloadbalancing:us-east-1:"$AWS_ACCOUNT_ID":targetgroup/genetraps-dnatoken-client/"$AWS_CLIENT_DNATOKEN_TARGET_GROUP --container-name client-dnatoken --container-port 8083 --ecs-profile genetraps

#./scripts/update-repo.sh

#dx rm -r -a /samples/

exit 0
