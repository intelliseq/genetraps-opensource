#!/bin/bash

INPUT_URL=$1
DESTINATION="/java-test"
PROJECT="project-F5qXjk008ZxZKP6XGYYzqbG7"

dx run url_fetcher \
            -i url="$INPUT_URL" \
            --destination "$PROJECT"":""$DESTINATION" \
            -y

