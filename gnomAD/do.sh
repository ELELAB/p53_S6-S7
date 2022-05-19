. /usr/local/envs/py37/bin/activate

./gnomad2csv TP53\
	-b GRCh37\
	-g gnomad_r2_1\
	-t 0.0001\
	-T exome\
	-o p53_gnomad.csv\
	--protein-only\
	--canonical-only

