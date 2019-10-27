#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

ctdna_frac = read_csv(file=url_gmaf_ctdna_frac, col_types = cols(.default = col_character()))  %>%
			 type_convert() %>%
			 group_by(gid, Tumor_Sample_Barcode) %>%
			 dplyr::summarize(ctdna_frac_pik3ca = max(ctdna_frac_pik3ca, na.rm=TRUE),
					   		  ctdna_frac_clonalmut = max(ctdna_frac_clonalmut, na.rm=TRUE),
					   		  ctdna_frac_long = max(ctdna_frac_long, na.rm=TRUE)) %>%
			 dplyr::rename(Guardant_ID = gid, MSK_ID = Tumor_Sample_Barcode)

facets_key = read_csv(file=url_facets_key, col_types = cols(.default = col_character()))  %>%
			 type_convert() %>%
			 dplyr::select(MSK_ID = `cBio Portal ID`) %>%
			 mutate(File_ID = MSK_ID)
			 
ctdna_frac = left_join(ctdna_frac, facets_key, by="MSK_ID")

cfdna_key = read_csv(file=url_cfdna_key, col_types = cols(.default = col_character()))  %>%
			type_convert() %>%
			mutate(Guardant_ID = paste0(Nov.ID, " ", sample_tx_status))

load(url_gmaf_to_annotate)
MAF = left_join(gmaf, cfdna_key, by="SampleId") %>%
	  left_join(cfdna_key, by="Guardant_ID") %>%
	  mutate(dbSNP_RS = dbSNP_RS.1, dbSNP_RS=ifelse(dbSNP_RS=="", "novel", dbSNP_RS),
	  		 SampleId = SampleId.x,
	  		 sample_tx_status = sample_tx_status.x,
	  		 Nov.ID = Nov.ID.x,
	  		 Cycle = Cycle.x,
	  		 Day = Day.x) %>%
	  dplyr::select(-Variant_Type.1,
	  				-HGVSp_Short.1,
	  		 		-dbSNP_RS.1,
	  		 		-Copy_number,
	  		 		-SampleId.x,
	  		 		-SampleId.y,
	  		 		-sample_tx_status.x,
	  		 		-sample_tx_status.y,
	  		 		-Nov.ID.x,
	  		 		-Nov.ID.y,
	  		 		-Cycle.x,
	  		 		-Cycle.y,
	  		 		-Day.x,
	  		 		-Day.y) %>%
	  left_join(ctdna_frac, by="Guardant_ID") %>%
	  filter(!is.na(File_ID))
	  
plot.0 = ggplot(MAF, aes(x = ctdna_frac_pik3ca, y = ctdna_frac_clonalmut, fill = Guardant_ID)) +
 		 geom_point(alpha = .8, size = 3, color = "black", shape = 21) +
 		 geom_abline(slope=1, intercept=0, color = "goldenrod3") + 
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
 		 labs(x="\nctDNA fraction PIK3CA\n", y="\nctDNA fraction clonal mutations\n") +
  		 coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
pdf(file=str_c(out_dir, "ctDNA_Fraction_PIK3CA_Clonal_Mutations.pdf"))
print(plot.0)
dev.off()

plot.0 = ggplot(MAF, aes(x = ctdna_frac_pik3ca, y = ctdna_frac_long, fill = Guardant_ID)) +
 		 geom_point(alpha = .8, size = 3, color = "black", shape = 21) +
 		 geom_abline(slope=1, intercept=0, color = "goldenrod3") + 
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
 		 labs(x="\nctDNA fraction PIK3CA\n", y="\nctDNA fraction `long`\n") +
  		 coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
pdf(file=str_c(out_dir, "ctDNA_Fraction_PIK3CA_Long.pdf"))
print(plot.0)
dev.off()

plot.0 = ggplot(MAF, aes(x = ctdna_frac_clonalmut, y = ctdna_frac_long, fill = Guardant_ID)) +
 		 geom_point(alpha = .8, size = 3, color = "black", shape = 21) +
 		 geom_abline(slope=1, intercept=0, color = "goldenrod3") + 
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
 		 labs(x="\nctDNA fraction clonal mutations\n", y="\nctDNA fraction `long`\n") +
  		 coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
pdf(file=str_c(out_dir, "ctDNA_Fraction_Clonal_Mutations_Long.pdf"))
print(plot.0)
dev.off()

oo = list()
sample_id = unique(MAF$Guardant_ID)
for (i in 1:length(sample_id)) {
	tmp = MAF %>%
		  filter(Guardant_ID == sample_id[i])
	alpha = unique(tmp$ctdna_frac_clonalmut)
	if (is.infinite(alpha)) {
		alpha = unique(tmp$ctdna_frac_long)
	}
	
	fhat = tmp$VAF/100
	n = tmp$mol_count
	tmp2 = read_tsv(file=paste0("../data/byl_facets_updated_042519/", tmp$File_ID[1], ".txt"), col_types = cols(.default = col_character()))  %>%
		   type_convert()
	Chr = tmp$Chromosome
	Pos = tmp$Start_Position
	index = NULL
	for (j in 1:length(Chr)) {
		indx = which(tmp2$chrom==Chr[j] & tmp2$loc.start<=Pos[j] & tmp2$loc.end>=Pos[j])
		if (length(indx)!=0) {
			index = c(index, indx)
		} else {
			index = c(index, NA)
		}
	}
	qt = tmp2$tcn.em[index]
	q1 = tmp2$lcn.em[index]
	q2 = qt-q1
	qt[is.na(qt)] = round(median(tmp2$tcn.em, na.rm=TRUE))
	q2[is.na(q2)] = ceiling(qt[is.na(q2)]/2)
	tmp3 = cancercellFraction(fhat, n, qt, q2, alpha, e=1e-5)
	tmp3 = cbind(tmp3, qt = qt, q2 = q2, fhat = fhat)
	tmp = cbind(tmp, tmp3)
	oo[[i]] = tmp
}
oo = do.call(rbind, oo)
write_tsv(oo, path=str_c(out_dir, "byl_maf_annotated_ccf.tsv"), append=FALSE, col_names=TRUE)

## compare CCF in ctDNA and tumor
xx = read_csv(file=url_gmaf_ctdna_frac_v1, col_types = cols(.default = col_character()))  %>%
	 type_convert() %>%
	 dplyr::select(-X1) %>%
	 filter(!is.na(gid)) %>%
	 mutate(Guardant_ID = gsub(pattern=" ", replacement="_", x=gid, fixed=TRUE)) %>%
	 mutate(UUID = paste0(Guardant_ID, ":", Hugo_Symbol, ":", HGVSp_Short))

yy = read_tsv(file=str_c(out_dir, "byl_maf_annotated_ccf.tsv"), col_types = cols(.default = col_character()))  %>%
	 type_convert() %>%
	 mutate(Guardant_ID = gsub(pattern=" ", replacement="_", x=Guardant_ID, fixed=TRUE)) %>%
	 mutate(UUID = paste0(Guardant_ID, ":", Hugo_Symbol, ":", HGVSp_Short))
	 
zz = full_join(xx %>% select(UUID, ccf, ctdna_frac_pik3ca, ctdna_frac_clonalmut), yy %>% select(UUID, cancer_cell_frac), by="UUID")

plot.0 = ggplot(zz, aes(x = ccf, y = cancer_cell_frac)) +
  		 geom_point(alpha = .8, size = 3, aes(fill = ctdna_frac_clonalmut, color = ctdna_frac_clonalmut), shape = 21) +
  		 scale_fill_gradient(low="yellow", high="blue") +
  		 scale_color_gradient(low="yellow", high="blue") +
  		 geom_abline(slope=1, intercept=0, color = "goldenrod3") + 
  		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
  		 labs(x="\nCCF (tumor)\n", y="\nCCF (cfDNA)\n") +
   		 coord_cartesian(xlim = c(0, 1), ylim=c(0, 1))
pdf(file=str_c(out_dir, "ctDNA_tumor_CCF.pdf"))
print(plot.0)
dev.off()

plot.0 = ggplot(oo, aes(x = VAF, y = cancer_cell_frac*100)) +
  		 geom_point(alpha = .5, aes(fill = ctdna_frac_clonalmut, color = ctdna_frac_clonalmut), shape=21, size=3) +
  		 geom_hline(yintercept=90, linetype=3, color="red") + 
  		 geom_vline(xintercept=10, linetype=3, color="red") + 
  		 scale_fill_gradient(low="yellow", high="blue") +
  		 scale_color_gradient(low="yellow", high="blue") +
  		 geom_abline(slope=1, intercept=0, color = "goldenrod3") + 
  		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
  		 labs(x="\nVAF (%)\n", y="\nCCF (%)\n") +
  		 scale_y_log10(
  		 ) +
  		 scale_x_log10(
  		 ) +
   		 coord_cartesian(xlim = c(0.01, 100), ylim=c(0.01, 100))
   		 
pdf(file=str_c(out_dir, "ctDNA_VAF_CCF.pdf"))
print(plot.0)
dev.off()

## compare CCF in ctDNA with Guardant's estimate
guardant_ccf = read_csv(file=url_guardant_ccf, col_types = cols(.default = col_character()))  %>%
 	  		   type_convert() %>%
 	  		   filter(!is.na(clon2)) %>%
 	  		   filter(type!="CNV") %>%
 	  		   mutate(uuid = paste0(patientrecordid, "_", gene, "_", alteration)) %>%
 	  		   mutate(clon2 = ifelse(clon2==0, 0.001, clon2)) %>%
 	  		   mutate(percentage = ifelse(percentage==0, 0.1, percentage))
 	  		   
plot.0 = ggplot(guardant_ccf, aes(x = percentage, y = clon2*100)) +
  		 geom_point(alpha = .85, shape=21, size=3, fill="salmon") +
  		 geom_abline(slope=1, intercept=0, color = "goldenrod3") + 
  		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
  		 labs(x="\nVAF (%)\n", y="\nGuardant estimated CCF (%)\n") +
  		 scale_y_log10(
  		 ) +
  		 scale_x_log10(
  		 ) +
   		 coord_cartesian(xlim = c(0.01, 100), ylim=c(0.01, 100))
   		 
pdf(file=str_c(out_dir, "ctDNA_VAF_CCF_Guardant.pdf"))
print(plot.0)
dev.off()

guardant_ccf = read_csv(file=url_guardant_ccf, col_types = cols(.default = col_character()))  %>%
 	  		   type_convert() %>%
 	  		   filter(!is.na(clon2)) %>%
 	  		   filter(type!="CNV") %>%
 	  		   mutate(uuid = paste0(patientrecordid, "_", gene, "_", alteration))
 	  		   
estimated_ccf = read_tsv(file=str_c(out_dir, "byl_maf_annotated_ccf.tsv"), col_types = cols(.default = col_character()))  %>%
				type_convert() %>%
				mutate(Mut_aa = ifelse(is.na(Mut_aa), "", Mut_aa)) %>%
				mutate(uuid = paste0(SampleId, "_", gene, "_", Mut_aa))

tmp = left_join(estimated_ccf, guardant_ccf, by="uuid") %>%
	  filter(!is.na(cancer_cell_frac)) %>%
	  filter(!is.na(clon2)) %>%
	  mutate(cancer_cell_frac = ifelse(cancer_cell_frac==0, 0.001, cancer_cell_frac)) %>%
	  mutate(clon2 = ifelse(clon2==0, 0.001, clon2)) %>%
	  mutate(gene_uuid = case_when(
	  			Hugo_Symbol=="PIK3CA" ~ "PIK3CA",
				Hugo_Symbol== "ESR1" ~ "ESR1",
				TRUE ~ "Other")) %>%
	 mutate(gene_uuid = factor(gene_uuid, levels=c("PIK3CA", "ESR1", "Other"), ordered=TRUE)) %>%
	 mutate(sample_uuid = factor(Nov.ID))
	 
fill = col_distinct(60)
names(fill) = str_pad(1:60, width=3, side="left", pad="0")
shape = c(21:23)
names(shape) = unique(tmp$gene_uuid)
				
plot.0 = ggplot(tmp, aes(x = cancer_cell_frac*100, y = clon2*100, shape = gene_uuid, fill = sample_uuid)) +
  		 geom_point(alpha = .8, size = 3) +
 		 scale_shape_manual(values = shape) +
 		 scale_fill_manual(values = fill) +
 		 geom_abline(slope=1, intercept=0, linetype=1, color="goldenrod3", size=1) +
 		 geom_smooth(method='lm', formula=y~x, linetype=1, se=TRUE, color="salmon", aes(x = cancer_cell_frac*100, y = clon2*100), inherit.aes = FALSE) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
 		 labs(x="CCF % (tumor/ctDNA)", y="CCF % (ctDNA only)") +
  		 scale_y_log10(
  		 	breaks = function(x) { c(0.1, 1, 10, 100) },
  			labels = function(x) { c("0.1", "1", "10", "100") }
  		 ) +
  		 scale_x_log10(
  		 	breaks = function(x) { c(0.1, 1, 10, 100) },
  			labels = function(x) { c("0.1", "1", "10", "100") }
  		 ) +
   		 coord_cartesian(xlim = c(0.09, 100), ylim=c(0.09, 100)) +
   		 guides(fill=guide_legend(title=c("Patient ID"), override.aes=list(shape=21)), shape=guide_legend(title=c("Gene Symbol")))
   		 
pdf(file=str_c(out_dir, "ctDNA_CCF_CCF_Guardant.pdf"), width=8)
print(plot.0)
dev.off()

export_x = tmp %>%
		   dplyr::select(ccf_tumor_ctdna = cancer_cell_frac,
		   				 ccf_ctdna_only = clon2,
		   				 patient_id = sample_uuid,
		   				 gene_symbol = gene_uuid)

write_tsv(export_x, path=str_c(out_dir, "ctDNA_CCF_CCF_Guardant.tsv"), append = FALSE, col_names = TRUE)

## compare CCF in ctDNA by Guardant and tumor
xx = read_csv(file=url_gmaf_ctdna_frac_v1, col_types = cols(.default = col_character()))  %>%
	 type_convert() %>%
	 select(-X1) %>%
	 filter(!is.na(gid)) %>%
	 mutate(Guardant_ID = gsub(pattern=" ", replacement="_", x=gid, fixed=TRUE)) %>%
	 mutate(UUID = paste0(Guardant_ID, ":", Hugo_Symbol, ":", HGVSp_Short))
	 
yy = tmp %>%
	 mutate(UUID = paste0(Nov.ID, "_", sample_tx_status, ":", Hugo_Symbol, ":", HGVSp_Short))
	 
tmp = left_join(yy, xx, by="UUID") %>%
	  filter(!is.na(ccf)) %>%
	  mutate(gene_uuid = case_when(
	  			Hugo_Symbol.x=="PIK3CA" ~ "PIK3CA",
				Hugo_Symbol.x== "ESR1" ~ "ESR1",
				TRUE ~ "Other")) %>%
	  mutate(gene_uuid = factor(gene_uuid, levels=c("PIK3CA", "ESR1", "Other"), ordered=TRUE)) %>%
	  mutate(sample_uuid = factor(Nov.ID))
	  
fill = col_distinct(60)
names(fill) = str_pad(1:60, width=3, side="left", pad="0")
shape = c(21:23)
names(shape) = unique(tmp$gene_uuid)
	
plot.0 = ggplot(tmp, aes(x = ccf*100, y = clon2*100, shape = gene_uuid, fill = sample_uuid)) +
  		 geom_point(alpha = .8, size = 3) +
 		 scale_shape_manual(values = shape) +
 		 scale_fill_manual(values = fill) +
 		 geom_abline(slope=1, intercept=0, linetype=1, color="goldenrod3", size=1) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
 		 labs(x="CCF % (tumor)", y="CCF % (ctDNA only)") +
   		 coord_cartesian(xlim = c(5, 100), ylim=c(5, 100)) +
   		 guides(fill=guide_legend(title=c("Patient ID"), override.aes=list(shape=21)), shape=guide_legend(title=c("Gene Symbol")))
   		 
pdf(file=str_c(out_dir, "ctDNA_CCF_Tumor_CCF_ctDNA_Guardant.pdf"), width=8)
print(plot.0)
dev.off()
