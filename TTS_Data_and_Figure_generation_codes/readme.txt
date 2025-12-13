In each of these folders containing recorded data related TTS/ITS has either or both of these files: "coarse_summary.csv", "coarse_instance_tts.csv".




'coarse_summary.csv' has the columns: 1)prefix 2)thld_delta 3)tau_star_val 4) median_TTS_val
 	
	1) prefix denotes the size of the problem under study: for Max3SAT "uf100", for 3R3X "3R3X_n112", etc.
	2) thld_delta denotes sampling period delta of the FN-annealer.
	3) tau_star_val denotes the computational runtime (Tcomp) chosen at that delta for finding instance wise TTS for all instances at that problem size.
	4) median_TTS_val denotes the median TTS across all instances at that problem size for given thld_delta and tau_star_val.

'coarse_instance_tts.csv' details the instance-wise TTS for each of the instances at each problem size for the given thld_delta and tau_star_val mentioned in 'coarse_summary.csv'. These instance wise TTS are used to calculate the median_TTS_val at that (thld_delta,tau_star_val) for that problem size. 'coarse_instance_tts.csv' has the following columns.
1) prefix 2)thld_delta 3)thld_max 4)instance 5)tau_star_val 6)succ_trials 7)total_trials 8)Ps_hat 9)TTS_at_tau_star

	1) prefix carries 2) thld_delta 5) tau_star_val carries the same meaning as mentioned for 'coarse_summary.csv'
	3) thld_max denotes the chosen value of 'C' for the FN-annealer parameter.
	4) instance denotes the problem instance at that problem size under study.
	7) total_trials denotes the total number of independent trials that instance is run for.
	6) success_trials denotes the number of hit/successes out of those total_trials
	8) Ps_hat denotes the success probability for that instance
	9) TTS_at_tau_star denotes the instance-wise TTS for that instance at the given parameters.


Each of these folders also contain the matlab code used for generating the TTS/ITS figures in the paper.

