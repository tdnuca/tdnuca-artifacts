#!/bin/bash

if [ "$#" -ne 1 ]
then
	echo "$0 [input]"
	exit 1
fi

INPUT=$1

sed -i "3s/^rsa/RSA/g" $INPUT
sed -i "5iFGX                        1             0" $INPUT
