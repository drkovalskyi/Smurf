#!/bin/bash

for DIR in `ls | grep ^RemoteGlidein`;
do
    echo crab -report -c $DIR
done

