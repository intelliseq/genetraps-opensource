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
ecs-cli configure profile --profile-name ECS_CLI_GENETRAPS_DEV --access-key $ECS_CLI_GENETRAPS_DEV_ACCESS_KEY --secret-key $ECS_CLI_GENETRAPS_DEV_KEY_ID
ecs-cli configure --cluster GENETRAPS_DEV --default-launch-type FARGATE --region us-east-1 --config-name ECS_CLI_GENETRAPS_DEV_CONF
./scripts/update-repo.sh

exit 0
