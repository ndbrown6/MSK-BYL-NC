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
	   filter(!(UID %in% nomut_nopik_x$sample_id)) %>%
	   filter(!duplicated(UUID))

GCCF_pre_tx = GCCF %>%
			  filter(sample_tx_status=="Pre" | (Nov.ID=="030" & sample_tx_status=="Pre") | (Nov.ID=="037" & sample_tx_status=="Pre")) %>%
			  dplyr::select(Hugo_Symbol, HGVSp_Short, cancer_cell_frac, fhat, Nov.ID)

GCCF_post_tx = GCCF %>%
			   filter(sample_tx_status=="Post" | (Nov.ID=="030" & sample_tx_status=="On") | (Nov.ID=="037" & sample_tx_status=="On")) %>%
			   dplyr::select(Hugo_Symbol, HGVSp_Short, cancer_cell_frac, fhat, Nov.ID)
			   
sample_id = unique(intersect(GCCF_pre_tx$Nov.ID, GCCF_post_tx$Nov.ID))

d1 = foreach (i=1:length(sample_id)) %dopar% {
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

d2 = foreach (i=1:length(sample_id)) %dopar% {
	gccf_pre = GCCF_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["fhat"]]
	names(gccf_pre) = GCCF_pre_tx %>%
					  filter(Nov.ID == sample_id[i]) %>%
					  mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					  .[["UUID"]]
	    
	gccf_post = GCCF_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["fhat"]]
	names(gccf_post) = GCCF_post_tx %>%
					   filter(Nov.ID == sample_id[i]) %>%
					   mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					   .[["UUID"]]
	
	gccf_x = gccf_y = vector(mode="numeric", length(unique(c(names(gccf_pre), names(gccf_post)))))
	names(gccf_x) = names(gccf_y) = unique(c(names(gccf_pre), names(gccf_post)))
	gccf_x[names(gccf_pre)] = gccf_pre
	gccf_y[names(gccf_post)] = gccf_post
	
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
	
	gccf_x2 = gccf_y2 = vector(mode="numeric", length(unique(c(names(gccf_pre), names(gccf_post)))))
	names(gccf_x2) = names(gccf_y2) = unique(c(names(gccf_pre), names(gccf_post)))
	gccf_x2[names(gccf_pre)] = gccf_pre
	gccf_y2[names(gccf_post)] = gccf_post
	
	gccf = data_frame(vaf_pre = gccf_x,
					  vaf_post = gccf_y,
					  ccf_pre = gccf_x2,
					  ccf_post = gccf_y2,
					  uuid = names(gccf_x)) %>%
					  mutate(vaf_pre  = ifelse(vaf_pre==0, 0.0001, vaf_pre)) %>%
					  mutate(vaf_post = ifelse(vaf_post==0, 0.0001, vaf_post)) %>%
					  mutate(ccf_pre  = ifelse(ccf_pre==0, 0.001, ccf_pre)) %>%
					  mutate(ccf_post = ifelse(ccf_post==0, 0.001, ccf_post))
					  
	fit = fun_lm_vaf(data=gccf)
	pr = predict(fit, gccf)
	d = log10(gccf$vaf_post/pr)
	names(d) = gccf$uuid
	return(invisible(d))
}

n = unlist(lapply(d1, length))
n = rep(sample_id, times=n)
delta = data_frame(delta_ccf = unlist(d1),
				   delta_vaf = unlist(d2),
				   sample_uuid = factor(n),
				   mutation_uuid = names(unlist(d1))) %>%
		mutate(mutation_uuid = unlist(lapply(strsplit(mutation_uuid, "_"), function(x) {x[1]}))) %>%
		mutate(mutation_uuid = case_when(
					mutation_uuid == "PIK3CA" ~ "PIK3CA",
					mutation_uuid == "ESR1" ~ "ESR1",
					TRUE ~ "Other")) %>%
		mutate(mutation_uuid = factor(mutation_uuid, levels=c("PIK3CA", "ESR1", "Other"), ordered=TRUE))
		
fill = col_distinct(60)
names(fill) = str_pad(1:60, width=3, side="left", pad="0")
shape = c(21:23)
names(shape) = unique(delta$mutation_uuid)

pdf(file=str_c(out_dir, "delta_CCF_VAF.pdf"), width=8)
plot.0 = ggplot(delta, aes(x = delta_ccf, y = delta_vaf, shape = mutation_uuid, fill = sample_uuid)) +
 		 geom_point(alpha = .8, size = 3) +
 		 scale_shape_manual(values = shape) +
 		 scale_fill_manual(values = fill) +
 		 geom_abline(slope=1, intercept=0, linetype=1, color="goldenrod3", size=1) +
 		 geom_smooth(method='lm', formula=y~x, linetype=1, se=TRUE, color="salmon", aes(x = delta_ccf, y = delta_vaf), inherit.aes = FALSE) +
 		 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15)) +
 		 labs(x=expression(Delta~CCF), y=expression(Delta~VAF)) +
   		 guides(fill=guide_legend(title=c("Patient ID"), override.aes=list(shape=21)), shape=guide_legend(title=c("Gene Symbol")))
print(plot.0)
dev.off()

export_x = delta %>%
		   dplyr::rename(patient_id = sample_uuid,
		   				 gene_symbol = mutation_uuid)
		   				 
write_tsv(export_x, path=str_c(out_dir, "delta_CCF_VAF.tsv"), append = FALSE, col_names = TRUE)

## partially collapsed mutations
fn = unique(unlist(lapply(d2, names)))
mz = matrix(NA, nrow=length(fn), ncol=length(d2))
rownames(mz) = fn
colnames(mz) = sample_id
for (i in 1:length(d2)) {
	fx = d2[[i]]
	mz[names(fx),i] = fx
}
gene_symbols = unlist(lapply(strsplit(rownames(mz), split="_"), function(x) { x[1] }))
uuid_gene_symbols = unique(gene_symbols)
uuid_2 = c("PIK3CA", "ARID1A")
uuid_3 = c("ESR1", "PTEN", "NF1")
uuid_gene_symbols = uuid_gene_symbols[!(uuid_gene_symbols %in% c(uuid_2, uuid_3))]

tmp.0 = NULL
for (i in 1:length(uuid_gene_symbols)) {
	index = which(gene_symbols==uuid_gene_symbols[i])
	if (length(index)==1) {
		tmp.0 = rbind(tmp.0, mz[index,,drop=FALSE])
	} else if (length(index)>1) {
		'row_min' <- function(x) {
			y = min(x, na.rm=TRUE)
			if (is.infinite(y)) {
				y = NA
			}
			return(y)
		}
		tmp.0 = rbind(tmp.0, apply(mz[index,,drop=FALSE], 2, row_min))
	}
}
rownames(tmp.0) = uuid_gene_symbols

tmp.1 = NULL
for (i in 1:length(uuid_2)) {
	index = which(gene_symbols==uuid_2[i])
	'row_min' <- function(x) {
		x = as.vector(na.omit(x))
		if (length(x)==0) {
			x = c(NA, NA)
		} else if (length(x)==1) {
			x = c(x, NA)
		} else {
			y = order(abs(x), decreasing=FALSE)
			x = c(x[y[1]], x[y[length(y)]])
		}
		return(x)
	}
	tmp.1 = rbind(tmp.1, apply(mz[index,,drop=FALSE], 2, row_min))
}
rownames(tmp.1) = paste0(rep(uuid_2, each=2), " (", rep(1:length(uuid_2), times=2), ")")

tmp.2 = NULL
for (i in 1:length(uuid_3)) {
	index = which(gene_symbols==uuid_3[i])
	'row_min' <- function(x) {
		x = as.vector(na.omit(x))
		if (length(x)==0) {
			x = c(NA, NA, NA)
		} else if (length(x)==1) {
			x = c(x, NA, NA)
		} else if (length(x)==2) {
			y = order(abs(x), decreasing=FALSE)
			x = c(x[y], NA)
		} else {
			y = order(abs(x), decreasing=FALSE)
			x = c(x[y[1:2]], x[y[length(y)]])
		}
		return(x)
	}
	tmp.2 = rbind(tmp.2, apply(mz[index,,drop=FALSE], 2, row_min))
}
rownames(tmp.2) = paste0(rep(uuid_3, each=3), " (", rep(1:length(uuid_3), times=3), ")")

mz = rbind(tmp.0, tmp.1, tmp.2)
colnames(mz) = sample_id

load(str_c(out_dir, "byl_ordered_margins.RData"))
mz = mz[gene_symbol, sample_id, drop=FALSE]

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed.pdf"), width=5, height=7)
corr_plot(corr=mz, method = "circle", type = "full", add = FALSE,
       	 col = col_fun(2), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
       	 outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
         addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
         na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=0.95)
dev.off()
