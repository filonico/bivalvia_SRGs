#!/bin/bash

# generate the list of fasta headers to be maintained one-by-one (format inconsistence)

# AIRC
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/A_irradians_concentricus/Airc_genomic_noIso.gff | sed -E 's/^.+ID=cds.evm.model./Airc_/; s/;Parent.+$//' | sort -u > 01b_assemblies_otherResources/A_irradians_concentricus/Airc_header_noIso.ls

# AMAR
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/A_marissinica/Amar_genomic_noIso.gff | sed -E 's/^.+ID=cds./Amar_/; s/;Parent.+$//' | sort -u > 01b_assemblies_otherResources/A_marissinica/Amar_header_noIso.ls

# APUR
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/A_purpuratus/Apur_genomic_noIso.gff | sed -E 's/^.+ID=cds.evm.model./Apur_/; s/;Parent.+$//; s/_/./2' | sort -u > 01b_assemblies_otherResources/A_purpuratus/Apur_header_noIso.ls

# CARI
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/C_ariakensis/Cari_genomic_noIso.gff | sed -E 's/^.+Parent=/Cari_/' | sort -u > 01b_assemblies_otherResources/C_ariakensis/Cari_header_noIso.ls

# CSIN
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/C_sinensis/Csin_genomic_noIso.gff | sed -E 's/^.+ID=cds.evm.model./Csin_/; s/;Parent.+$//; s/_/./g; s/\./_/' | sort -u > 01b_assemblies_otherResources/C_sinensis/Csin_header_noIso.ls

# HBIA
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/H_bialata/Hbia_genomic_noIso.gff | sed -E 's/^.+Parent=Upi/Hbia_/' | sort -u > 01b_assemblies_otherResources/H_bialata/Hbia_header_noIso.ls

# MMAR
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/M_margaritifera/Mmar_genomic_noIso.gff | sed -E 's/^.+Parent=/Mmar_/' | sort -u > 01b_assemblies_otherResources/M_margaritifera/Mmar_header_noIso.ls

# MNER
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/M_nervosa/Mner_genomic_noIso.gff | sed -E 's/^.+Parent=/Mner_/' | sort -u > 01b_assemblies_otherResources/M_nervosa/Mner_header_noIso.ls

# MPHI
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/M_philippinarum/Mphi_genomic_noIso.gff | sed -E 's/^.+Mph_/Mphi_/; s/-mRNA.+$//; s/_/./2; s/-/./g' | sort -u > 01b_assemblies_otherResources/M_philippinarum/Mphi_header_noIso.ls

# PVIR
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/P_viridis/Pvir_genomic_noIso.gff | sed -E 's/^.+Parent=pvir_/Pvir_/; s/\..+$//' | sort -u > 01b_assemblies_otherResources/P_viridis/Pvir_header_noIso.ls

# SBRO
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/S_broughtonii/Sbro_genomic_noIso.gff | sed -E 's/^.+Parent=/Sbro_/' | sort -u > 01b_assemblies_otherResources/S_broughtonii/Sbro_header_noIso.ls

# SCON
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/S_constricta/Scon_genomic_noIso.gff | sed -E 's/^.+evm.model./Scon_/' | sort -u > 01b_assemblies_otherResources/S_constricta/Scon_header_noIso.ls

# SGLO
grep -P '\t'CDS'\t' 01b_assemblies_otherResources/S_glomerata/Sglo_genomic_noIso.gff | sed -E 's/^.+ID=/Sglo_/; s/-.+$//' | sort -u > 01b_assemblies_otherResources/S_glomerata/Sglo_header_noIso.ls

for i in 01b_assemblies_otherResources/*/*noIso*gff; do
	
	DIR=$(dirname $i) &&
	spID=$(basename $i | awk -F "_" '{print $1}') &&
	OUT_LIST=$(echo "$DIR"/"$spID"_header_noIso.ls) &&
	
	# extract survived sequences from CDS fastas	
	python3 scripts/08_extract_sequences_from_fasta.py -l $OUT_LIST -f "$DIR"/*rn.fna -o "$DIR"/"$spID"_cds_noIso.fna

	# extract survived sequences from protein fastas
        python3 scripts/08_extract_sequences_from_fasta.py -l $OUT_LIST -f "$DIR"/*rn.faa -o "$DIR"/"$spID"_protein_noIso.faa


done
