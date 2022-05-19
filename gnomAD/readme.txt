#use HGNC notation for the gene name 
#modify do.sh script with gene name
#have gnomad2csv in the folder where you run  
./do.sh
cat p53_gnomad.csv | awk -F "\"*,\"*" '{print $11}' > p53_gnomad_mut_list.txt 

