#!/bin/bash

for r in {1..4}
do
	sbatch NMA_error_single_dep.sh $r
	sbatch NMA_error_prev_dep.sh $r
done

