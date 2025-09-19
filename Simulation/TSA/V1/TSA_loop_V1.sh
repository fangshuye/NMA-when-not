#!/bin/bash

for r in {1..3}
do
	sbatch TSA_V1.sh $r
done

