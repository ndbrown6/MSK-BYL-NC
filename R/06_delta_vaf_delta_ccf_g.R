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
			   filter(!is.na(Nov.ID))

pre_ids = unique(cfdna_key$Nov.ID[which(cfdna_key$sample_tx_status == 'Pre' & !cfdna_key$SampleId %in% c(duplicate_id$sample_id, nomut_nopik_x$sample_id))])
post_ids = unique(cfdna_key$Nov.ID[which(cfdna_key$sample_tx_status == 'Post' & !cfdna_key$SampleId %in% c(duplicate_id$sample_id, nomut_nopik_x$sample_id))])

pre_ids = c(pre_ids, "030", "037")
post_ids = c(post_ids, "030", "037")

sample_id = unique(intersect(pre_ids, post_ids))
sample_id = sample_id[sample_id != "042"]

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
	names(gmaf_pre) = gmaf_pre_tx %>% filter(Nov.ID == sample_id[i]) %>% mutate(UUID = paste0(gene, "_", alteration)) %>% .[["UUID"]]
	    
	gmaf_post = gmaf_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["clon2"]]
	names(gmaf_post) = gmaf_post_tx %>% filter(Nov.ID == sample_id[i]) %>% mutate(UUID = paste0(gene, "_", alteration)) %>% .[["UUID"]]
	
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

cfdna_key = read_csv(file=url_cfdna_key, col_types = cols(.default = col_character()))  %>%
	   		type_convert() %>%
	   		mutate(Guardant_ID = paste0(Nov.ID, " ", sample_tx_status))
			
load(url_gmaf_to_annotate)
oncokb_maf = read_tsv(file=url_gmaf_oncokb_annotated, col_types = cols(.default = col_character()))  %>%
	   		 type_convert() %>%
			 mutate(oncokb_effect = ifelse(oncogenic=="", 0, 1)) %>%
			 mutate(oncokb_effect = ifelse(is.na(oncokb_effect), 0, 1)) %>%
			 select(oncokb_effect)
gmaf = cbind(gmaf, oncokb_maf)

gmaf_pre_tx = gmaf %>%
			  filter(sample_tx_status=="Pre" | (Nov.ID=="030" & sample_tx_status=="Pre") | (Nov.ID=="037" & sample_tx_status=="Pre")) %>%
			  filter(Nov.ID %in% sample_id) %>%
			  dplyr::select(Hugo_Symbol, HGVSp_Short, VAF, Nov.ID) %>%
			  mutate(VAF = VAF/100)

gmaf_post_tx = gmaf %>%
			   filter(sample_tx_status=="Post" | (Nov.ID=="030" & sample_tx_status=="On") | (Nov.ID=="037" & sample_tx_status=="On")) %>%
			   filter(Nov.ID %in% sample_id) %>%
			   dplyr::select(Hugo_Symbol, HGVSp_Short, VAF, Nov.ID) %>%
			   mutate(VAF = VAF/100)
			   
d2 = foreach (i=1:length(sample_id)) %dopar% {
	gmaf_pre = gmaf_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["VAF"]]
	names(gmaf_pre) = gmaf_pre_tx %>% filter(Nov.ID == sample_id[i]) %>% mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>% .[["UUID"]]
	    
	gmaf_post = gmaf_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["VAF"]]
	names(gmaf_post) = gmaf_post_tx %>% filter(Nov.ID == sample_id[i]) %>% mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>% .[["UUID"]]
	
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
					TRUE ~ "Other")) %>%
		mutate(gene_uuid = factor(gene_uuid, levels=c("PIK3CA", "ESR1", "Other"), ordered=TRUE)) %>%
		mutate(sample_uuid = factor(sample_uuid))
		
fill = col_distinct(60)
names(fill) = str_pad(1:60, width=3, side="left", pad="0")
shape = c(21:23)
names(shape) = c("Other", "ESR1", "PIK3CA")

pdf(file=str_c(out_dir, "delta_CCF_VAF_by_Guardant.pdf"), width=8)
plot.0 = ggplot(delta, aes(x = x, y = y, shape = gene_uuid, fill = sample_uuid)) +
 		 geom_point(alpha = .8, size = 3) +
 		 scale_shape_manual(values = shape) +
 		 scale_fill_manual(values = fill) +
 		 geom_abline(slope=1, intercept=0, linetype=1, color="goldenrod3", size=1) +
 		 geom_smooth(method='lm', formula=y~x, linetype=1, se=TRUE, color="salmon", aes(x = x, y = y), inherit.aes = FALSE) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
 		 labs(x=expression(Delta~"CCF (ctDNA only)"), y=expression(Delta~"VAF")) +
 		 guides(fill=guide_legend(title=c("Patient ID"), override.aes=list(shape=21)), shape=guide_legend(title=c("Gene Symbol")))
print(plot.0)
dev.off()
