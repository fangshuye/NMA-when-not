#!/bin/bash

for r in {1..4}
do
	sbatch NMA_error_single.sh $r
done

