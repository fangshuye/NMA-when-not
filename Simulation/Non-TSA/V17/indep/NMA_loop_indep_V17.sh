#!/bin/bash

for r in {1..4}
do
	sbatch NMA_error_single_indep.sh $r
	sbatch NMA_error_prev_indep.sh $r
done

