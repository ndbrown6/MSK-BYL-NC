#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

GCCF = read_tsv(file=str_c(out_dir, "byl_maf_annotated_ccf.tsv"), col_types = cols(.default = col_character()))  %>%
	   type_convert() %>%
	   mutate(UID = paste0(Nov.ID, " ", Cycle)) %>%
	   mutate(UUID = paste0(SampleId, "_", Hugo_Symbol, ":", Chromosome, ":", Start_Position, ":", Reference_Allele, ">", Tumor_Seq_Allele2)) %>%
	   filter(!(UID %in% nomut_nopik_y$sample_id)) %>%
	   filter(!duplicated(UUID))

GCCF_pre_tx = GCCF %>%
			  filter(sample_tx_status=="Pre" | (Nov.ID=="030" & sample_tx_status=="Pre") | (Nov.ID=="037" & sample_tx_status=="Pre")) %>%
			  dplyr::select(Hugo_Symbol, HGVSp_Short, cancer_cell_frac, fhat, Nov.ID)

GCCF_post_tx = GCCF %>%
			   filter(sample_tx_status=="Post" | (Nov.ID=="030" & sample_tx_status=="On") | (Nov.ID=="037" & sample_tx_status=="On")) %>%
			   dplyr::select(Hugo_Symbol, HGVSp_Short, cancer_cell_frac, fhat, Nov.ID)

sample_id = unique(intersect(GCCF_pre_tx$Nov.ID, GCCF_post_tx$Nov.ID))			   
d2 = foreach (i=1:length(sample_id)) %dopar% {
	gccf_pre = GCCF_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["cancer_cell_frac"]]
	names(gccf_pre) = GCCF_pre_tx %>%
					  filter(Nov.ID == sample_id[i]) %>%
					  mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					  .[["UUID"]]
	    
	gccf_post = GCCF_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["cancer_cell_frac"]]
	names(gccf_post) = GCCF_post_tx %>%
					   filter(Nov.ID == sample_id[i]) %>%
					   mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					   .[["UUID"]]
	
	gccf_x = gccf_y = vector(mode="numeric", length(unique(c(names(gccf_pre), names(gccf_post)))))
	names(gccf_x) = names(gccf_y) = unique(c(names(gccf_pre), names(gccf_post)))
	gccf_x[names(gccf_pre)] = gccf_pre
	gccf_y[names(gccf_post)] = gccf_post
	
	gccf = data_frame(pre_treatment = gccf_x,
					  post_treatment = gccf_y,
					  uuid = names(gccf_x)) %>%
					  mutate(pre_treatment = ifelse(pre_treatment==0, 0.001, pre_treatment)) %>%
					  mutate(post_treatment = ifelse(post_treatment==0, 0.001, post_treatment))
					  
	fit = fun_lm_ccf(data=gccf %>% rename(ccf_pre=pre_treatment, ccf_post=post_treatment))
	pr = predict(fit, gccf %>% rename(ccf_pre=pre_treatment, ccf_post=post_treatment))
	d = log10(gccf$post_treatment/pr)
	names(d) = gccf$uuid
	return(invisible(d))
}

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
			   filter(!is.na(Nov.ID))

gmaf_pre_tx = guardant_ccf %>%
			  filter(sample_tx_status=="Pre" | (Nov.ID=="030" & sample_tx_status=="Pre") | (Nov.ID=="037" & sample_tx_status=="Pre")) %>%
			  filter(Nov.ID %in% sample_id) %>%
			  dplyr::select(gene, alteration, clon2, Nov.ID)

gmaf_post_tx = guardant_ccf %>%
			   filter(sample_tx_status=="Post" | (Nov.ID=="030" & sample_tx_status=="On") | (Nov.ID=="037" & sample_tx_status=="On")) %>%
			   filter(Nov.ID %in% sample_id) %>%
			   dplyr::select(gene, alteration, clon2, Nov.ID)
			   

d1 = foreach (i=1:length(sample_id)) %dopar% {
 	gmaf_pre = gmaf_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["clon2"]]
	names(gmaf_pre) = gmaf_pre_tx %>%
					  filter(Nov.ID == sample_id[i]) %>%
					  mutate(UUID = paste0(gene, "_", alteration)) %>%
					  .[["UUID"]]
	    
	gmaf_post = gmaf_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["clon2"]]
	names(gmaf_post) = gmaf_post_tx %>%
					   filter(Nov.ID == sample_id[i]) %>%
					   mutate(UUID = paste0(gene, "_", alteration)) %>%
					   .[["UUID"]]
	
	gmaf_x = gmaf_y = vector(mode="numeric", length(unique(c(names(gmaf_pre), names(gmaf_post)))))
	names(gmaf_x) = names(gmaf_y) = unique(c(names(gmaf_pre), names(gmaf_post)))
	gmaf_x[names(gmaf_pre)] = gmaf_pre
	gmaf_y[names(gmaf_post)] = gmaf_post
 	
 	gccf = data_frame(vaf_pre = gmaf_x,
 					  vaf_post = gmaf_y,
 					  uuid = names(gmaf_x)) %>%
 					  mutate(vaf_pre  = ifelse(vaf_pre==0, 1e-5, vaf_pre)) %>%
 					  mutate(vaf_post = ifelse(vaf_post==0, 1e-5, vaf_post))
 					  
 	fit = fun_lm_vaf_2(data=gccf, sample_id = sample_id[i])
 	pr = predict(fit, gccf)
 	d = log10(gccf$vaf_post/pr)
 	d = ifelse(abs(d)<1e-10, 0, d)
 	names(d) = gccf$uuid
 	return(invisible(d))
}

d3 = foreach (i=1:length(sample_id)) %dopar% {
	featureNames = intersect(names(d1[[i]]), names(d2[[i]]))
	x = d1[[i]][featureNames]
	y = d2[[i]][featureNames]
	z = data_frame(x=x,
				   y=y,
				   sample_uuid = sample_id[i],
				   gene_uuid = unlist(lapply(strsplit(featureNames, "_", fixed=TRUE), function(x) {x[1]})))
	return(z)
}

delta = do.call(rbind, d3) %>%
		mutate(gene_uuid = case_when(
					gene_uuid == "PIK3CA" ~ "PIK3CA",
					gene_uuid == "ESR1" ~ "ESR1",
					gene_uuid == "PTEN" ~ "PTEN",
					gene_uuid == "NF1" ~ "NF1",
					gene_uuid == "TP53" ~ "TP53",
					TRUE ~ "Other")) %>%
		mutate(gene_uuid = factor(gene_uuid)) %>%
		mutate(sample_uuid = factor(sample_uuid))
		
fill = col_distinct(26)
names(fill) = as.vector(unique(delta$sample_uuid))
shape = c(21:25, 8)
names(shape) = unique(delta$gene_uuid)

pdf(file=str_c(out_dir, "delta_CCF_CCF_by_Guardant.pdf"))
plot.0 = ggplot(delta, aes(x = x, y = y, shape = gene_uuid, fill = sample_uuid)) +
 		 geom_point(alpha = .8, size = 3) +
 		 scale_shape_manual(values = shape) +
 		 scale_fill_manual(values = fill) +
 		 geom_abline(slope=1, intercept=0, linetype=1, color="goldenrod3", size=1) +
 		 geom_smooth(method='lm', formula=y~x, linetype=1, se=TRUE, color="salmon", aes(x = x, y = y), inherit.aes = FALSE) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
 		 labs(x=expression(Delta~"CCF (Guardant)"), y=expression(Delta~"CCF (Absolute)")) +
 		 coord_cartesian(xlim = c(-4,4), ylim = c(-4,4)) +
 		 guides(fill=guide_legend(title=c("Patient ID"), override.aes=list(shape=21)), shape=guide_legend(title=c("Gene Symbol")))
print(plot.0)
dev.off()
