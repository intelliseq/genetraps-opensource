#!/bin/bash

JOB_ID="$1"

dx describe "$JOB_ID" | grep "State\|Failure message"


