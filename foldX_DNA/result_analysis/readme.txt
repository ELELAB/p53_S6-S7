#aim: create table s2. 
#requirements:

	ddG_overview.csv (from results)
	saturation_mutagenesis_overview.csv
	selected_neutral_list.txt

ln -s ../compare_saturation_mutagenesis/saturation_mutagenesis_overview.csv
ln -s ../selected_neutral/selected_neutral_list.txt

module load python/3.7/modulefile
python3 compare_ddGstability_ddGDNA.py
