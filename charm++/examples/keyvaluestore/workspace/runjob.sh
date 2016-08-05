#!/bin/sh
qsub -A ATPESC2016 -t 5 -n 4 --mode script script.vesta
