#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

guardant_ccf = read_csv(file=url_guardant_ccf, col_types = cols(.default = col_character()))  %>%
 	  		   type_convert() %>%
 	  		   filter(!is.na(clon2)) %>%
 	  		   filter(type!="CNV") %>%
 	  		   rename(SampleId = patientrecordid)

cfdna_key = read_csv(file=url_cfdna_key, col_types = cols(.default = col_character()))  %>%
	   		type_convert() %>%
	   		mutate(Guardant_ID = paste0(Nov.ID, " ", sample_tx_status))

guardant_ccf = left_join(guardant_ccf, cfdna_key, by="SampleId") %>%
			   filter(!is.na(SampleId)) %>%
			   filter(!is.na(Nov.ID)) %>%
			   mutate(UUID = paste0(Nov.ID, "_", Cycle, "_", gene, "_", alteration))


load(url_gmaf_to_annotate)
oncokb_maf = read_tsv(file=url_gmaf_oncokb_annotated, col_types = cols(.default = col_character()))  %>%
 	  		 type_convert() %>%
			 mutate(oncokb_effect = ifelse(oncogenic=="", 0, 1)) %>%
			 mutate(oncokb_effect = ifelse(is.na(oncokb_effect), 0, 1)) %>%
			 select(oncokb_effect)
gmaf = bind_cols(gmaf, oncokb_maf) %>%
	   mutate(UUID = paste0(Nov.ID, "_", Cycle, "_", gene, "_", Mut_aa))


guardant_ccf = left_join(guardant_ccf, gmaf, by="UUID") %>%
			   mutate(UUID = paste0(Nov.ID.x, "_", gene.x, "_", HGVSp_Short)) %>%
			   filter(UUID %in% paste0(max_vaf$Sample_ID, "_", max_vaf$HGVSp_Short))

rm(list=ls(all=TRUE))
source("config.R")
gmaf_ctdna_frac = read_csv(file=url_gmaf_ctdna_frac_v1, col_types = cols(.default = col_character()))  %>%
				  type_convert() %>%
				  dplyr::select(-X1) %>%
				  filter(Hugo_Symbol=="PIK3CA" & is.na(clonal_status=="clonal") %>%
				  filter(!is.na(gid)) %>%
				  mutate(Guardant_ID = gsub(pattern=" ", replacement="_", x=gid, fixed=TRUE)) %>%
				  mutate(UUID = paste0(Guardant_ID, ":", Hugo_Symbol, ":", HGVSp_Short))

rm(list=ls(all=TRUE))
source("config.R")
load(url_gmaf_to_annotate)
export_x = gmaf %>%
		   filter(Nov.ID == "009" & sample_tx_status == "Pre") %>%
		   dplyr::select(gene_symbol = Hugo_Symbol,
		   				 chromosome = Chromosome,
		   				 start_position = Start_Position,
		   				 end_position = End_Position,
		   				 reference_allele = Reference_Allele,
		   				 alternate_allele = Tumor_Seq_Allele2,
		   				 context_3nt = flanking_bps)
		   				 
mutation_summary = export_x %>%
				   mutate(patient_id = "009_Prex_Tx",
				   		  Type = "SNV") %>%
				   dplyr::select(Sample = patient_id,
				   				 CHROM = chromosome,
				   				 POS = start_position,
				   				 REF = reference_allele,
				   				 ALT = alternate_allele,
				   				 Type)
vcf = preprocessInput_snv(input_data = mutation_summary,
                          ensgene = ensgene,
                          reference_genome = ref_genome)
                          
plot_96_spectrum(vcf, sample.col = "Sample",  file = str_c(out_dir, "mutation_signatures_009_Pre_Tx.pdf"))
write_tsv(export_x, path=str_c(out_dir, "mutation_signatures_009_Pre_Tx.tsv"), append = FALSE, col_names = TRUE)

export_x = gmaf %>%
		   filter(Nov.ID == "040" & sample_tx_status == "Pre") %>%
		   filter(Reference_Allele!="-") %>%
		   filter(Tumor_Seq_Allele2!="-") %>%
		   dplyr::select(gene_symbol = Hugo_Symbol,
		   				 chromosome = Chromosome,
		   				 start_position = Start_Position,
		   				 end_position = End_Position,
		   				 reference_allele = Reference_Allele,
		   				 alternate_allele = Tumor_Seq_Allele2,
		   				 context_3nt = flanking_bps)
		   				 
mutation_summary = export_x %>%
				   mutate(patient_id = "040_Prex_Tx",
				   		  Type = "SNV") %>%
				   dplyr::select(Sample = patient_id,
				   				 CHROM = chromosome,
				   				 POS = start_position,
				   				 REF = reference_allele,
				   				 ALT = alternate_allele,
				   				 Type)
vcf = preprocessInput_snv(input_data = mutation_summary,
                          ensgene = ensgene,
                          reference_genome = ref_genome)
                          
plot_96_spectrum(vcf, sample.col = "Sample",  file = str_c(out_dir, "mutation_signatures_040_Pre_Tx.pdf"))
write_tsv(export_x, path=str_c(out_dir, "mutation_signatures_040_Pre_Tx.tsv"), append = FALSE, col_names = TRUE)
