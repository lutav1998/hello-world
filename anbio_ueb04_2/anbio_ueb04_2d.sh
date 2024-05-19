#!/bin/bash

# Reis-Transkriptom "all.cdna" herunterladen und benennen
wget -O http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.cdna all.cdna;

# Indexierung namens "o_sativa" mit all.cdna durchführen 
kallisto index -i o_sativa.idx all.cdna;

# Kallisto Loop (Quantifizierung)
for i in fastq/files/path/*.fastq.gz; do	
	# Ordner nicht berücksichtigen						
	if [ -f $i]; then	
		# Quantifizierung aller FastQ Dateien entsprechend den .gz Dateien
		kallisto quant -i o_sativa.idx -o . --single -l 200 -s 40 $i
	fi
done;

 