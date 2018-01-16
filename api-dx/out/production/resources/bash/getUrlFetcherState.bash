#!/bin/bash

JOB_ID="$1"
PROJECT="project-F5qXjk008ZxZKP6XGYYzqbG7"

dx describe \
            --project "$PROJECT" \
            "$JOB_ID" | grep "State\|Failure message"


