#!/bin/bash

for instance in ../instances/*; do
	echo $instance: >> target.txt
done
