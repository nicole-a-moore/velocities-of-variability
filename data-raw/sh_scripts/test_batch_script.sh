#!/bin/bash
for chosen_one in {1..21}
do

sbatch test_job_script.sh $chosen_one

done