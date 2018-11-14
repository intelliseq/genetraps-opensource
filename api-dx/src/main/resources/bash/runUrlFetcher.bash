#!/bin/bash

INPUT_URL=$1
SAMPLEID=$2
DESTINATION="samples/$SAMPLEID/rawdata"

dx run url_fetcher \
            -i url="$INPUT_URL" \
            --folder="$DESTINATION" \
            -y

