#!/bin/bash

# rename gff files for each species
zcat 01b_assemblies_otherResources/S_glomerata/S.glomerata_genemodels.gff3.gz > 01b_assemblies_otherResources/S_glomerata/Sglo_genomic.gff

zcat 01b_assemblies_otherResources/A_irradians_concentricus/Argopecten_irradians_concentricus.gff.gz > 01b_assemblies_otherResources/A_irradians_concentricus/Airc_genomic.gff

zcat 01b_assemblies_otherResources/A_marissinica/Ama_v1.0.gff.gz > 01b_assemblies_otherResources/A_marissinica/Amar_genomic.gff

zcat 01b_assemblies_otherResources/A_purpuratus/Scallop.gff.gz > 01b_assemblies_otherResources/A_purpuratus/Apur_genomic.gff

zcat 01b_assemblies_otherResources/C_ariakensis/Chr_genome_final_gene.gff3.gz | sed -E '/CDS.+Parent=EVM0015005.1/ s/  11493781/\t  11493781/g' | sed '/exon/ s/  11512101/\t  11512101/' | sed -E 's/11515623 11515849/11515623\t  11515849/' > 01b_assemblies_otherResources/C_ariakensis/Cari_genomic.gff

zcat 01b_assemblies_otherResources/C_sinensis/Cyclina_sinensis_genome.gff3.gz > 01b_assemblies_otherResources/C_sinensis/Csin_genomic.gff

zcat 01b_assemblies_otherResources/H_bialata/Upi_annotation_v4.gff3.gz > 01b_assemblies_otherResources/H_bialata/Hbia_genomic.gff

zcat 01b_assemblies_otherResources/M_margaritifera/Mma_Braker2_Agat_Id_Ovl_StartStopFixed_CDS100.gff3.gz > 01b_assemblies_otherResources/M_margaritifera/Mmar_genomic.gff

zcat 01b_assemblies_otherResources/M_nervosa/ReferenceGenomeSequences/AllSupport.gff3.gz > 01b_assemblies_otherResources/M_nervosa/Mner_genomic.gff

zcat 01b_assemblies_otherResources/M_philippinarum/Modiolus_philippinarum.gff.gz > 01b_assemblies_otherResources/M_philippinarum/Mphi_genomic.gff

zcat 01b_assemblies_otherResources/P_viridis/P.viridis_genome-annotation.gff3.gz > 01b_assemblies_otherResources/P_viridis/Pvir_genomic.gff

zcat 01b_assemblies_otherResources/S_broughtonii/EVM.final.gene.gff3.gz > 01b_assemblies_otherResources/S_broughtonii/Sbro_genomic.gff

zcat 01b_assemblies_otherResources/S_constricta/Sco_annotation.gff3.gz > 01b_assemblies_otherResources/S_constricta/Scon_genomic.gff

zcat 01b_assemblies_otherResources/S_glomerata/S.glomerata_genemodels.gff3.gz > 01b_assemblies_otherResources/S_glomerata/Sglo_genomic.gff

# remove isoforms with agat
for i in 01b_assemblies_otherResources/*/*_genomic.gff; do

        agat_sp_keep_longest_isoform.pl -gff $i -o "${i::-4}"_noIso.gff &&
        rm -r $i

done
