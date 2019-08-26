FROM ubuntu:18.04

RUN echo "cd /tmp" > cromwell.sh
RUN echo "wget http://resouces.intelliseq.pl/public/intelliseqngs/workflows/cromwell/cromwell.jar" >> cromwell.sh
RUN echo "wget http://resouces.intelliseq.pl/public/intelliseqngs/workflows/cromwell/aws.cfg" >> cromwell.sh
RUN echo 'sed -i -e "s+BUCKET_TAG+$BUCKET_TAG+" -e "s+QUEUE_TAG+$QUEUE_TAG+" aws.cfg' >> cromwell.sh
RUN echo "echo starting cromwell server" >> cromwell.sh
RUN echo "java -Dconfig.file=aws.cfg -jar cromwell.jar server" >> cromwell.sh

ENTRYPOINT ["sh", "./cromwell.sh"]
