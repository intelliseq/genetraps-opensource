LOG_PREFIX="===== TRAVIS LOG ===== "

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
echo $LOG_PREFIX"api-dx: building..."
gradle build docker -p api-dx/
echo $LOG_PREFIX"api-dx: running docker..."
docker run -d -p 8080:8080 -e "DNANEXUS_TOKEN="$DNANEXUS_TOKEN -t pl.intelliseq.genetraps.api.dx/api-dx:latest
echo $LOG_PREFIX"api-dx: waiting for service..."
./scripts/wait-for-service.sh localhost:8080/hello 60
echo $LOG_PREFIX"api-dx: checking for error..."
check
echo $LOG_PREFIX"api-dx: testing..."
echo $(curl localhost:8080/touch)

### BUILDING CLINET_EXPLORARE ###
echo $LOG_PREFIX"client-explorare: compiling..."
npm --prefix ./client-explorare install ./client-explorare
echo $LOG_PREFIX"clien-explorare: testing..."
npm --prefix client-explorare run test
echo $LOG_PREFIX"client-explorare: dockerizing..."
docker build -t pl.intelliseq.genetraps.client.explorare/client-explorare:latest client-explorare
echo $LOG_PREFIX"client-explorare: running docker..."
docker run -d -p 8081:80 pl.intelliseq.genetraps.client.explorare/client-explorare:latest
echo $LOG_PREFIX"client-explorare: waiting for service..."
./scripts/wait-for-service.sh localhost:8081/hello 60
echo $LOG_PREFIX"client-explorare: checking for error..."
check
exit 0
