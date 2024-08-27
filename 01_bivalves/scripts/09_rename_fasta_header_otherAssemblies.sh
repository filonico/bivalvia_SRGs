#!/bin/bash

# AIRC
zcat 01b_assemblies_otherResources/A_irradians_concentricus/Argopecten_irradians_concentricus_pep.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/evm\.model\./Airc_/' > 01b_assemblies_otherResources/A_irradians_concentricus/Airc_protein_rn.faa
zcat 01b_assemblies_otherResources/A_irradians_concentricus/Argopecten_irradians_concentricus_cds.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/evm\.model\./Airc_/' > 01b_assemblies_otherResources/A_irradians_concentricus/Airc_cds_rn.fna

# AMAR
zcat 01b_assemblies_otherResources/A_marissinica/Ama_v1.0_pep.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>/>Amar_/' > 01b_assemblies_otherResources/A_marissinica/Amar_protein_rn.faa
zcat 01b_assemblies_otherResources/A_marissinica/Ama_v1.0_transcript.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>/>Amar_/' > 01b_assemblies_otherResources/A_marissinica/Amar_cds_rn.fna

# APUR
zcat 01b_assemblies_otherResources/A_purpuratus/Scallop.pep.fasta.gz |  bash scripts/05_onelinerizer.sh - | sed -E 's/evm.model./Apur_/; s/_/./2' > 01b_assemblies_otherResources/A_purpuratus/Apur_protein_rn.faa
zcat 01b_assemblies_otherResources/A_purpuratus/Scallop.cds.fasta.gz |  bash scripts/05_onelinerizer.sh - | sed -E 's/evm.model./Apur_/; s/_/./2' > 01b_assemblies_otherResources/A_purpuratus/Apur_cds_rn.fna

# CARI
zcat 01b_assemblies_otherResources/C_ariakensis/Chr_genome_final_gene.gff3.pep.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>/>Cari_/' > 01b_assemblies_otherResources/C_ariakensis/Cari_protein_rn.faa
zcat 01b_assemblies_otherResources/C_ariakensis/Chr_genome_final_gene.gff3.cds.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>/>Cari_/' > 01b_assemblies_otherResources/C_ariakensis/Cari_cds_rn.fna

# CSIN
zcat 01b_assemblies_otherResources/C_sinensis/Cyclina_sinensis_genome.pep.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/evm\.model\.//; s/_/./g; s/^>/>Csin_/' > 01b_assemblies_otherResources/C_sinensis/Csin_protein_rn.faa
zcat 01b_assemblies_otherResources/C_sinensis/Cyclina_sinensis_genome.cds.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/evm\.model\.//; s/_/./g; s/^>/>Csin_/' > 01b_assemblies_otherResources/C_sinensis/Csin_cds_rn.fna

# HBIA
zcat 01b_assemblies_otherResources/H_bialata/Upi_proteins_v4.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>Upi/>Hbia_/' > 01b_assemblies_otherResources/H_bialata/Hbia_protein_rn.faa
zcat 01b_assemblies_otherResources/H_bialata/Upi_cds_v4.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>Upi/>Hbia_/' > 01b_assemblies_otherResources/H_bialata/Hbia_cds_rn.fna

# MMAR
zcat 01b_assemblies_otherResources/M_margaritifera/Mma_Braker2_agat_final_proteins.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>/>Mmar_/' > 01b_assemblies_otherResources/M_margaritifera/Mmar_protein_rn.faa
zcat 01b_assemblies_otherResources/M_margaritifera/Mma_Braker2_agat_final_transcripts.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>/>Mmar_/' > 01b_assemblies_otherResources/M_margaritifera/Mmar_cds_rn.fna

# MNER
zcat 01b_assemblies_otherResources/M_nervosa/ReferenceGenomeSequences/AllSupport.prot.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+//; s/^>/>Mner_/' > 01b_assemblies_otherResources/M_nervosa/Mner_protein_rn.faa
zcat 01b_assemblies_otherResources/M_nervosa/ReferenceGenomeSequences/AllSupport.cds.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+//; s/^>/>Mner_/' > 01b_assemblies_otherResources/M_nervosa/Mner_cds_rn.fna

# MPHI
zcat 01b_assemblies_otherResources/M_philippinarum/Modiolus_philippinarum.pep.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>Mph/>Mphi/; s/_/./2; s/-/./' > 01b_assemblies_otherResources/M_philippinarum/Mphi_protein_rn.faa
zcat 01b_assemblies_otherResources/M_philippinarum/Modiolus_philippinarum.cds.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>Mph/>Mphi/; s/_/./2; s/-/./' > 01b_assemblies_otherResources/M_philippinarum/Mphi_cds_rn.fna

# PVIR
zcat 01b_assemblies_otherResources/P_viridis/P.viridis_genemodel.pep.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>p/>P/; s/\.p1//' > 01b_assemblies_otherResources/P_viridis/Pvir_protein_rn.faa
zcat 01b_assemblies_otherResources/P_viridis/P.viridis_genemodel.mrna.fa.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>p/>P/' > 01b_assemblies_otherResources/P_viridis/Pvir_cds_rn.fna

# SBRO
zcat 01b_assemblies_otherResources/S_broughtonii/EVM.final.gene.gff3.pep.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>/>Sbro_/' > 01b_assemblies_otherResources/S_broughtonii/Sbro_protein_rn.faa
zcat 01b_assemblies_otherResources/S_broughtonii/EVM.final.gene.gff3.cds.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/ .+$//; s/^>/>Sbro_/' > 01b_assemblies_otherResources/S_broughtonii/Sbro_cds_rn.fna

# SCON
# CDSs ARE MISSING
zcat 01b_assemblies_otherResources/S_constricta/Sco_proteins.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/evm.model./Scon_/; s/_/./2' > 01b_assemblies_otherResources/S_constricta/Scon_protein.rn.faa

# SGLO
zcat 01b_assemblies_otherResources/S_glomerata/S.glomerata_genemodels.pep.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>/>Sglo_/; s/-.+$//' > 01b_assemblies_otherResources/S_glomerata/Sglo_protein_rn.faa
zcat 01b_assemblies_otherResources/S_glomerata/S.glomerata_genemodels.trans.fasta.gz | bash scripts/05_onelinerizer.sh - | sed -E 's/^>/>Sglo_/; s/-.+$//' > 01b_assemblies_otherResources/S_glomerata/Sglo_cds_rn.fna

