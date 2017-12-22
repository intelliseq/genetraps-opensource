#echo $AWS_ACCOUNT_ID
#!/bin/bash
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
ecs-cli push $TAG_API_DX
ecs-cli push $TAG_CLIENT_EXPLORARE
#echo $TAG_CLIENT_EXPLORARE 
cat aws-conf/docker-compose-template.yml | sed 's@clientExplorareImageTag@'"$TAG_CLIENT_EXPLORARE"'@' > docker-compose.yml
ecs-cli compose --project-name genetraps-client-explorare -f docker-compose.yml --ecs-params ./aws-conf/ecs-params.yml service up --target-group-arn "arn:aws:elasticloadbalancing:us-east-1:"$AWS_ACCOUNT_ID":targetgroup/genetraps-explorare-client/"$AWS_GENETRAPS_TARGET_GROUP --container-name client-explorare --container-port 8081
