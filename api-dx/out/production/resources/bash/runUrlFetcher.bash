#!/bin/bash

INPUT_URL=$1
SAMPLE_NUMBER=$2
DESTINATION="samples/$SAMPLE_NUMBER/rawdata"

dx run url_fetcher \
            -i url="$INPUT_URL" \
            --folder="$DESTINATION" \
            -y

