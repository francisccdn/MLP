#!/bin/bash

echo "--TSP Benchmark--"

make

k=1
for instance in instances/*; do
	echo $instance >> ./benchmark/bm.txt

	echo "Running $instance"
	echo "Instance $k of 67" 

	for i in {1..10}; do
		./tsp ${instance} | grep 'COST\|TIME' | awk "{print $1}" >> ./benchmark/bm.txt
	done

	k=$(($k + 1))
done

echo "-" >> ./benchmark/bm.txt

echo "Running bm.py to compute averages..."

cd benchmark
python3 bm.py

echo "Finishing up summary..."
python3 summarycount.py
cd ..

echo "Benchmark completed."