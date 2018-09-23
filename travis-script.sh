LOG_PREFIX="===== TRAVIS LOG ===== "
LOG_APP=" no-app: "
#
echo $LOG_PREFIX "PYTHON VERSION: " `/usr/bin/python --version`
echo $LOG_PREFIX "JAVA VERSION: " `javac -version`
echo $LOG_PREFIX"Starting travis script"
if [ -z "$HELLO" ]; then source scripts/get-secrets.sh; fi
echo $LOG_PREFIX"Checking first five letters of token: "`echo $DNANEXUS_TOKEN_TEST | cut -c1-5`
dx login --token $DNANEXUS_TOKEN_TEST --noprojects
dx select genetraps-test

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

### CONFIGURING ECS ###
echo $LOG_PREFIX $LOG_APP "ecs-cli configuring"
ecs-cli configure --cluster genetraps --default-launch-type FARGATE --region us-east-1 --config-name genetraps

#############################
### BUILDING CLIENT_INDEX ###
#############################

LOG_APP="client-index: "
echo $LOG_PREFIX $LOG_APP "setting tag"
CLIENT_INDEX_CHECKSUM=`find client-index -type f -exec md5sum {} \; | sort -k 2 | md5sum | sed 's/  -//g'`
CLIENT_INDEX_TAG=$AWS_ACCOUNT_ID".dkr.ecr."$AWS_REGION".amazonaws.com/genetraps-client-index:"$CLIENT_INDEX_CHECKSUM
echo $LOG_PREFIX $LOG_APP $CLIENT_INDEX_TAG
echo $LOG_PREFIX $LOG_APP "checking if image exists"
CLIENT_INDEX_EXISTS=`aws ecr list-images --repository-name genetraps-client-index --profile genetraps | grep $CLIENT_INDEX_CHECKSUM | wc -l`
echo $LOG_PREFIX $LOG_APP "exists:" $EXISTS
#docker pull $CLIENT_INDEX_TAG
if [ $CLIENT_INDEX_EXISTS -eq 1 ]; then
    echo $LOG_PREFIX"docker image already exists"
else
    echo $LOG_PREFIX"building new image"
    docker build client-index/ -t $CLIENT_INDEX_TAG -q
    echo $LOG_PREFIX $LOG_APP "client tag: " $CLIENT_INDEX_TAG
    docker push $CLIENT_INDEX_TAG | cat
    cat aws-conf/docker-compose-template.yml | \
        sed 's@appTag@'"client-index"'@' | \
        sed 's@imageTag@'"$CLIENT_INDEX_TAG"'@' | \
        sed 's@portTag@'"$ECS_CLI_PORT_CLIENT_INDEX"'@g' | \
        sed 's@prefixTag@client-index-log@' > docker-compose.yml
    echo $LOG_PREFIX $LOG_APP "ecs-cli composing client-index"
    ecs-cli compose \
        --project-name genetraps-client-index \
        -f docker-compose.yml \
        --ecs-params ./aws-conf/ecs-params.yml \
        service up \
            --target-group-arn "arn:aws:elasticloadbalancing:"$AWS_REGION":"$AWS_ACCOUNT_ID":targetgroup/client-index-target-group/"$ECS_CLI_TG_CLIENT_INDEX \
            --container-name client-index \
            --container-port $ECS_CLI_PORT_CLIENT_INDEX \
            --aws-profile genetraps \
            --timeout 1
fi

#############################
### BUILDING API_SECURITY ###
#############################

LOG_APP="api-security: "
echo $LOG_PREFIX $LOG_APP "setting tag"
API_SECURITY_CHECKSUM=`find api-security -type f -exec md5sum {} \; | sort -k 2 | grep -P "/src/|Readme" | md5sum | sed 's/  -//g'`
API_SECURITY_TAG=$AWS_ACCOUNT_ID".dkr.ecr."$AWS_REGION".amazonaws.com/genetraps-api-security:"$API_SECURITY_CHECKSUM
echo $LOG_PREFIX $LOG_APP $API_SECURITY_TAG
#docker pull $API_SECURITY_TAG
API_SECURITY_EXISTS=`aws ecr list-images --repository-name genetraps-api-security --profile genetraps | grep $API_SECURITY_CHECKSUM | wc -l`
if [ $API_SECURITY_EXISTS -eq 1 ]; then
    echo $LOG_PREFIX"docker image already exists"
else
    echo $LOG_PREFIX"building new image"
    gradle build -q -p api-security 2>&1 | grep -v "WARNING:"
    cp `ls api-security/build/libs/api-security*` api-security/build/libs/app.jar
    docker build api-security/ -t $API_SECURITY_TAG -q
    echo $LOG_PREFIX $LOG_APP "security tag: " $API_SECURITY_TAG
    docker push $API_SECURITY_TAG | cat
    cat aws-conf/docker-compose-template.yml | \
        sed 's@appTag@'"api-security"'@' | \
        sed 's@imageTag@'"$API_SECURITY_TAG"'@' | \
        sed 's@portTag@'"$ECS_CLI_PORT_API_SECURITY"'@g' | \
        sed 's@prefixTag@api-security-log@' > docker-compose.yml
    echo "    environment:" >> docker-compose.yml
    echo "      - AURORA_GENETRAPS_CLIENT_PASSWD" >> docker-compose.yml
    echo "      - AURORA_GENETRAPS_CLIENT_LOGIN" >> docker-compose.yml
    echo $LOG_PREFIX $LOG_APP "ecs-cli composing api-security"
    ecs-cli compose \
        --project-name genetraps-api-security \
        -f docker-compose.yml \
        --ecs-params ./aws-conf/ecs-params.yml \
        service up \
            --target-group-arn "arn:aws:elasticloadbalancing:"$AWS_REGION":"$AWS_ACCOUNT_ID":targetgroup/api-security-target-group/"$ECS_CLI_TG_API_SECURITY \
            --container-name api-security \
            --container-port $ECS_CLI_PORT_API_SECURITY \
            --aws-profile genetraps \
            --timeout 1
fi

#############################
### BUILDING API_DX #########
#############################

LOG_APP="api-dx: "
echo $LOG_PREFIX $LOG_APP "setting tag"
API_DX_CHECKSUM=`find api-dx -type f -exec md5sum {} \; | sort -k 2 | grep -P "/src/|Readme" | md5sum | sed 's/  -//g'`
API_DX_TAG=$AWS_ACCOUNT_ID".dkr.ecr."$AWS_REGION".amazonaws.com/genetraps-api-dx:"$API_DX_CHECKSUM
echo $LOG_PREFIX $LOG_APP $API_DX_TAG
#docker pull $API_DX_TAG
API_DX_EXISTS=`aws ecr list-images --repository-name genetraps-api-dx --profile genetraps | grep $API_DX_CHECKSUM | wc -l`
if [ $API_DX_EXISTS -eq 1 ]; then
    echo $LOG_PREFIX"docker image already exists"
else
    echo $LOG_PREFIX"building new image"
    gradle build -q -p api-dx
    cp `ls api-dx/build/libs/api-dx*` api-dx/build/libs/app.jar
    docker build api-dx/ -t $API_DX_TAG -q
    echo $LOG_PREFIX $LOG_APP "dx tag: " $API_DX_TAG
    docker push $API_DX_TAG | cat
    cat aws-conf/docker-compose-template.yml | \
        sed 's@appTag@'"api-dx"'@' | \
        sed 's@imageTag@'"$API_DX_TAG"'@' | \
        sed 's@portTag@'"$ECS_CLI_PORT_API_DX"'@g' | \
        sed 's@prefixTag@api-security-log@' > docker-compose.yml
    echo "    environment:" >> docker-compose.yml
    echo "      - AURORA_GENETRAPS_CLIENT_PASSWD" >> docker-compose.yml
    echo "      - AURORA_GENETRAPS_CLIENT_LOGIN" >> docker-compose.yml
    echo "      - DNANEXUS_TOKEN_TEST" >> docker-compose.yml
    echo $LOG_PREFIX $LOG_APP "ecs-cli composing api-dx"
    ecs-cli compose \
        --project-name genetraps-api-dx \
        -f docker-compose.yml \
        --ecs-params ./aws-conf/ecs-params.yml \
        service up \
            --target-group-arn "arn:aws:elasticloadbalancing:"$AWS_REGION":"$AWS_ACCOUNT_ID":targetgroup/api-dx-target-group/"$ECS_CLI_TG_API_DX \
            --container-name api-dx \
            --container-port $ECS_CLI_PORT_API_DX \
            --aws-profile genetraps \
            --timeout 1
fi

################################
### BUILDING CLIENT_DX #########
################################

LOG_APP="client-dx: "
echo $LOG_PREFIX $LOG_APP "setting tag"
CLIENT_DX_CHECKSUM=`find client-dx -type f -exec md5sum {} \; | sort -k 2 | md5sum | sed 's/  -//g'`
CLIENT_DX_TAG=$AWS_ACCOUNT_ID".dkr.ecr."$AWS_REGION".amazonaws.com/genetraps-client-dx:"$CLIENT_DX_CHECKSUM
echo $LOG_PREFIX $LOG_APP $CLIENT_DX_TAG
#docker pull $CLIENT_DX_TAG
CLIENT_DX_EXISTS=`aws ecr list-images --repository-name genetraps-client-dx --profile genetraps | grep $CLIENT_DX_CHECKSUM | wc -l`
if [ $CLIENT_DX_EXISTS -eq 1 ]; then
    echo $LOG_PREFIX"docker image already exists"
else
    echo $LOG_PREFIX"building new image"
    docker build client-dx/ -t $CLIENT_DX_TAG -q
    echo $LOG_PREFIX $LOG_APP "dx tag: " $CLIENT_DX_TAG
    docker push $CLIENT_DX_TAG | cat
    cat aws-conf/docker-compose-template.yml | \
        sed 's@appTag@'"client-dx"'@' | \
        sed 's@imageTag@'"$CLIENT_DX_TAG"'@' | \
        sed 's@portTag@'"$ECS_CLI_PORT_CLIENT_DX"'@g' | \
        sed 's@prefixTag@api-security-log@' > docker-compose.yml
    echo $LOG_PREFIX $LOG_APP "ecs-cli composing client-dx"
    ecs-cli compose \
        --project-name genetraps-client-dx \
        -f docker-compose.yml \
        --ecs-params ./aws-conf/ecs-params.yml \
        service up \
            --target-group-arn "arn:aws:elasticloadbalancing:"$AWS_REGION":"$AWS_ACCOUNT_ID":targetgroup/client-dx-target-group/"$ECS_CLI_TG_CLIENT_DX \
            --container-name client-dx \
            --container-port $ECS_CLI_PORT_CLIENT_DX \
            --aws-profile genetraps \
            --timeout 1
fi




### BUILDING API_EXPLORARE
### BUILDING API_EXPLORARE ###
#APP="api-explorare"
#PORT=$ECS_CLI_PORT_API_EXPLORARE
#TG=$ECS_CLI_TG_API_EXPLORARE
#LOG_APP=$APP": "
#echo $LOG_PREFIX $LOG_APP "building..."
#gradle build -p $APP/
#echo $LOG_PREFIX $LOG_APP "setting tag"
#CHECKSUM=`find $APP -type f -exec md5sum {} \; | sort -k 2 | md5sum | sed 's/  -//g'`
#TAG=""
#TAG=$AWS_ACCOUNT_ID".dkr.ecr."$AWS_REGION".amazonaws.com/genetraps/"$APP":"$CHECKSUM
#echo $LOG_PREFIX $LOG_APP $TAG
#docker build $APP"/" -t $TAG -q
#echo $LOG_PREFIX $LOG_APP " tag=" $TAG
#docker run -d -p $PORT:$PORT -t $TAG
#echo $LOG_PREFIX $LOG_APP "waiting for service..."
#./scripts/wait-for-service.sh localhost:$PORT/hello 60
#echo $LOG_PREFIX $LOG_APP "checking for error..."
#check
#docker push $TAG
#cat aws-conf/docker-compose-template.yml | sed 's@appTag@'"$APP"'@' | sed 's@imageTag@'"$TAG"'@' | sed 's@portTag@'"$PORT"'@g' | sed 's@prefixTag@'"$APP"'-log@' > docker-compose.yml
#echo $LOG_PREFIX $LOG_APP "ecs-cli composing "$APP
#ecs-cli compose --project-name genetraps-$APP -f docker-compose.yml --ecs-params ./aws-conf/ecs-params.yml service up --target-group-arn "arn:aws:elasticloadbalancing:"$AWS_REGION":"$AWS_ACCOUNT_ID":targetgroup/$APP-target-group/"$TG --container-name $APP --container-port $PORT --aws-profile genetraps

### BUILDING CLIENT_DX ###
#LOG_APP="client-dx: "
#echo $LOG_PREFIX $LOG_APP "building..."
#gradle build docker -p client-dx/
#echo $LOG_PREFIX $LOG_APP "running docker..."
#docker run -d -p 8086:8086 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN_TEST -t pl.intelliseq.genetraps.api.dx/client-dx:latest
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

#TAG_CLIENT_DX=`ls -alR --full-time client-dx/ -Ibin -Ibuild -I.* | sha1sum | cut -d " " -f1`
#TAG_CLIENT_EXPLORARE=`ls -alR --full-time client-explorare/ -Inode_modules -Idist -Ietc -I.* | tr -s ' ' | cut -f1-5,9 | sha1sum | cut -d " " -f1`
#TAG_CLIENT_DNATOKEN=`ls -alR --full-time client-dnatoken/ -Inode_modules -Idist -Ietc -I.* | tr -s ' ' | cut -f1-5,9 | sha1sum | cut -d " " -f1`
#echo $TAG_CLIENT_DX
#echo $TAG_CLIENT_EXPLORARE
#IMAGE_CLIENT_DX=pl.intelliseq.genetraps.api.dx/client-dx:latest
#IMAGE_CLIENT_EXPLORARE=pl.intelliseq.genetraps.client.explorare/client-explorare:latest
#IMAGE_CLIENT_DNATOKEN=pl.intelliseq.genetraps.client.dnatoken/client-dnatoken:latest
#echo $IMAGE_CLIENT_EXPLORARE
#TAG_CLIENT_DX=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-client-dx:"$TAG_CLIENT_DX
#TAG_CLIENT_EXPLORARE=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-client-explorare:"$TAG_CLIENT_EXPLORARE
#TAG_CLIENT_DNATOKEN=$AWS_ACCOUNT_ID".dkr.ecr.us-east-1.amazonaws.com/genetraps-client-dnatoken:"$TAG_CLIENT_DNATOKEN
#docker tag $IMAGE_CLIENT_DX $TAG_CLIENT_DX
#docker tag $IMAGE_CLIENT_EXPLORARE $TAG_CLIENT_EXPLORARE
#docker tag $IMAGE_CLIENT_DNATOKEN $TAG_CLIENT_DNATOKEN
#ecs-cli push $TAG_CLIENT_DX --ecs-profile genetraps
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
