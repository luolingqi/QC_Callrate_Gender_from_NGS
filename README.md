# QC_Callrate_Gender_from_NGS
# Check sample call-rate and gender from BAM-Readcount using 96 Biomarks SNPs

This app is used as an internal QC, by predicting sample call rate and gender using only the 96 Biomarks SNPs. It reads in both the sample sheet in the run folder and the output from "BAM-Readcount", which extracts the position specific count information. The call rate is defined as the percentage of the 96 SNP sites which pass the minimum reads cutoff (default: 5). The rule of gender calling is exactly same as in the Biomarks User Guide. This app generates an output file to report the call rate and gender (PASS/FAIL), as well as a log file recording read depth at each single SNP site.
