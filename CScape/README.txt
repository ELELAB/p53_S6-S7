
Contact: Freja Dahl Hede, s174333@student.dtu.dk

This folder contains scripts and results for mutation driver status for p53 project. 

FOLDERS:
### scripts ###
convert_cancermuts_data.R		: Converts cancermuts metatable to get hgvs format on separate lines.
					  Usage: ./convert_aggregation_data.R
					  Output: /results/cancermuts_converted_hg19.csv

### results ###
cancermuts_converted_hg19.csv		: Cancermuts metatable converted from hgvs to separate columns. 
cancermuts_converted_hg19_CScape.csv	: Same as file above, but only with chr, pos, ref, mut and without header:
					  cut -f1,2,3,4 -d "," cancermuts_converted_hg19.csv | tail -n+2 > cancermuts_converted_hg19_CScape.csv
CScape_results.csv			: Output from CScape when run on cancermuts_converted_hg19_CScape.csv. Command:
					  cscape_somatic_query CScape/results/cancermuts_converted_hg19_CScape.csv -o /data/user/shared_projects/p53_jmb_2021/CScape/results/CScape_results.csv
cancermuts_CScape_MCAP_Revel.csv	: Cancermuts data (for the residues with SNVs) with CScape, M-CAP and Revel scores.
cancermutsDBD_CScape_MCAP_Revel.csv	: Same as above, but only in DBD (residue 91-289)

