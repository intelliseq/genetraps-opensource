#!/bin/bash

SERVICE_URL="$1"
TIMEOUT="$2"

echo "Waiting for endpoint $1"

start=`date +%s`
while [ -z ${DISCOVERY_SERVICE_READY} ]; do
  echo "Waiting for endpoint $1"
  if [ "$(curl --silent $SERVICE_URL 2>&1 | grep -q '\"status\":\"up\"'; echo $?)" = 0 ]; then
      DISCOVERY_SERVICE_READY=true;
  fi
  end=`date +%s`
  #echo $((end-start))
  if [ $((end-start)) -gt $TIMEOUT ]; then
      exit 1
  fi
  sleep 2
done

echo "Service started"
