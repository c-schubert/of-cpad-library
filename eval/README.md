## Evaluation of results
Evaluation script for cpa files written by the cpad library. 

Script |  Function
------------ | -------------
eval_seriell_cpas | Read and filter (get rid of some numerical errors in phase areas etc.) cpa results files from serial ran post processing utility or solver and write result file
eval_par_cpas | Read and filter multiple cpa files (for all used processor) and consoildate to results file
eval_cpas_final | Read the result file and do things with it ...