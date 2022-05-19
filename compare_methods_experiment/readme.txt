cp ../thermomutDB/tableExport_DBD_single.csv .

mutatex

cd md_mutatex

ln -s ../mutatex/zinc_bound/2XWRa_91-289/md/replicate1/CHARMM22star/saturation_rename/results/mutation_ddgs/final_averages/ final_averages

cd ..

cd xray_mutatex

ln -s ../mutatex/zinc_bound/2XWRa_91-289/xray/saturation/results/mutation_ddgs/final_averages/ final_averages

cd ..

cd pdbredo_mutatex

ln -s ../mutatex/zinc_bound/2XWRa_91-289/pdbredo/saturation/results/mutation_ddgs/final_averages/ final_averages

Rosetta:

cd xray_rosetta
cp ../rosetta_analysis/zinc_bound/2XWRa_91-289/exp/ref2015/aggregate/ddg_mutations_aggregate.csv .

cd ..

cd pdbredo_rosetta
cp ../rosetta_analysis/zinc_bound/2XWRa_91-289/pdbredo/ref2015/aggregate/ddg_mutations_aggregate.csv .

cd ..

#notice that the pdbredo_rosetta was updated 03-Mar-2021 after a new aggregation

module load python/3.7/modulefile

python3 scripts/compare_performance.py

#output 
#ddG_compare.csv

#plotting

python3 scripts/experimental_overview.py

#outputs two scatterplots describing the experimental data:
#	avg_themomut_ddg.pdf
#	differencies_in_thermomuts.pdf

python3 scripts/plot_compare.py

#outputs two plots to desribe the comparison between the methods and experiment
#	regplot.pdf
#	heatmap.pdf
