### login
echo $DOCKERHUB_PASSWORD | docker login -u $DOCKERHUB_USERNAME --password-stdin
### build docker
cat fastqc/Dockerfile | docker build - -t intelliseq/fastqc
### test
docker run -v "$PWD"/../../test-data/fastq/:/tmp/data -v /tmp/:/tmp/out intelliseq/fastqc -o /tmp/out /tmp/data/capn3.1.fq.gz

