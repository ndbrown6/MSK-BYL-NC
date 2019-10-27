#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

cfdna_key = read_csv(file=url_cfdna_key, col_types = cols(.default = col_character()))  %>%
	   		type_convert() %>%
	   		mutate(Guardant_ID = paste0(Nov.ID, " ", sample_tx_status))
			
load(url_gmaf_to_annotate)
oncokb_maf = read_tsv(file=url_gmaf_oncokb_annotated, col_types = cols(.default = col_character()))  %>%
	   		 type_convert() %>%
			 mutate(oncokb_effect = ifelse(oncogenic=="", 0, 1)) %>%
			 select(oncokb_effect)
gmaf = bind_cols(gmaf, oncokb_maf)

pre_ids = unique(cfdna_key$Nov.ID[which(cfdna_key$sample_tx_status == 'Pre' & !cfdna_key$SampleId %in% c(duplicate_id$sample_id,nomut_nopik_x$sample_id))])
post_ids = unique(cfdna_key$Nov.ID[which(cfdna_key$sample_tx_status == 'Post' & !cfdna_key$SampleId %in% c(duplicate_id$sample_id,nomut_nopik_x$sample_id))])

pre_ids = c(pre_ids, "030", "037")
post_ids = c(post_ids, "030", "037")

sample_id = unique(intersect(pre_ids, post_ids))
sample_id = sample_id[sample_id != "042"]
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

pdf(file=str_c(out_dir, "ctDNA_VAF_Pre_by_Post.pdf"), width=6, height=6)
pb <- txtProgressBar(min=0, max=length(sample_id), style=3)
for (i in 1:length(sample_id)) {
	setTxtProgressBar(pb, i)
	
	gmaf_pre = gmaf_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["VAF"]]
	names(gmaf_pre) = gmaf_pre_tx %>%
					  filter(Nov.ID == sample_id[i]) %>%
					  mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					  .[["UUID"]]
	    
	gmaf_post = gmaf_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["VAF"]]
	names(gmaf_post) = gmaf_post_tx %>%
					   filter(Nov.ID == sample_id[i]) %>%
					   mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					   .[["UUID"]]
	
	gmaf_x = gmaf_y = vector(mode="numeric", length(unique(c(names(gmaf_pre), names(gmaf_post)))))
	names(gmaf_x) = names(gmaf_y) = unique(c(names(gmaf_pre), names(gmaf_post)))
	gmaf_x[names(gmaf_pre)] = gmaf_pre
	gmaf_y[names(gmaf_post)] = gmaf_post
	
	gccf = data_frame(pre_treatment = gmaf_x,
					  post_treatment = gmaf_y,
					  uuid = names(gmaf_x)) %>%
					  mutate(pre_treatment = ifelse(pre_treatment==0, 1e-5, pre_treatment)) %>%
					  mutate(post_treatment = ifelse(post_treatment==0, 1e-5, post_treatment)) %>%
					  mutate(sample_uuid = sample_id[i])
    plot.0 = ggplot(gccf, aes(x = pre_treatment*100, y = post_treatment*100, label=uuid)) +
 		     geom_point(alpha = .8, size = 3, color = "black", shape = 21, fill = "salmon") +
 		 	 geom_abline(slope=1, intercept=0, color = "goldenrod3") +
 		 	 geom_text_repel() +
 		 	 scale_x_log10(
 		 		breaks = function(x) { c(1e-3, 0.1, 1, 10, 100) },
  				labels = function(x) { c("1E-3", "0.1", "1", "10", "100") }
 		 	 ) +
 		 	 scale_y_log10(
 		 		breaks = function(x) { c(1e-3, 0.1, 1, 10, 100) },
  				labels = function(x) { c("1E-3", "0.1", "1", "10", "100") }
 		 	 ) +
 		 	 annotation_logticks() +
 		 	 theme_bw(base_size = 17) + 
 		 	 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
 		 	 labs(x="\nVAF (% pre-treatment)\n", y="\nVAF (% post-treatment)\n") +
  		 	 coord_cartesian(xlim = 100*c(1e-5, 1), ylim = 100*c(1e-5, 1)) +
  		 	 facet_wrap(~sample_uuid)
	print(plot.0)
}
close(pb)
dev.off()

pdf(file=str_c(out_dir, "ctDNA_VAF_Pre_by_Post_L.pdf"), width=5, height=6)
pb <- txtProgressBar(min=0, max=length(sample_id), style=3)
for (i in 1:length(sample_id)) {
	setTxtProgressBar(pb, i)
	
	gccf_pre = gmaf_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["VAF"]]
	names(gccf_pre) = gmaf_pre_tx %>%
					  filter(Nov.ID == sample_id[i]) %>%
					  mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					  .[["UUID"]]
	    
	gccf_post = gmaf_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["VAF"]]
	names(gccf_post) = gmaf_post_tx %>%
					   filter(Nov.ID == sample_id[i]) %>%
					   mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					   .[["UUID"]]
	
	gccf_x = gccf_y = vector(mode="numeric", length(unique(c(names(gccf_pre), names(gccf_post)))))
	names(gccf_x) = names(gccf_y) = unique(c(names(gccf_pre), names(gccf_post)))
	gccf_x[names(gccf_pre)] = gccf_pre
	gccf_y[names(gccf_post)] = gccf_post
	
	gccf = data_frame(vaf = c(gccf_x, gccf_y),
					  time = factor(c(rep("Pre", length(gccf_x)),rep("Post", length(gccf_y))), levels=c("Pre", "Post")),
					  uuid = rep(names(gccf_x), 2)) %>%
					  mutate(vaf = ifelse(vaf==0, 0.001, vaf)) %>%
					  mutate(sample_uuid = sample_id[i])
    plot.0 = ggplot(gccf, aes(x = time, y = vaf*100, group=uuid)) +
    		 geom_line() +
    		 geom_point() +
 		 	 scale_y_log10(
  		 		breaks = function(x) { c(0.1, 1, 10, 100) },
   				labels = function(x) { c("0.1", "1", "10", "100") }
  		 	 ) +
 		 	 annotation_logticks(side="l") +
 		 	 theme_bw(base_size = 17) + 
 		 	 theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), legend.position = "none") +
 		 	 labs(x=" ", y="\nVAF (%)\n") +
  		 	 coord_cartesian(xlim = c(1.45,1.55), ylim = 100*c(0.001, 1)) +
  		 	 facet_wrap(~sample_uuid)
	print(plot.0)
}
dev.off()

d = foreach (i=1:length(sample_id)) %dopar% {
	gmaf_pre = gmaf_pre_tx %>% 
	    	   filter(Nov.ID == sample_id[i]) %>%
	    	   .[["VAF"]]
	names(gmaf_pre) = gmaf_pre_tx %>%
					  filter(Nov.ID == sample_id[i]) %>%
					  mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
					  .[["UUID"]]
	    
	gmaf_post = gmaf_post_tx %>% 
	    		filter(Nov.ID == sample_id[i]) %>%
	    		.[["VAF"]]
	names(gmaf_post) = gmaf_post_tx %>%
					   filter(Nov.ID == sample_id[i]) %>%
					   mutate(UUID = paste0(Hugo_Symbol, "_", HGVSp_Short)) %>%
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
 
## fully collapsed mutations
fn = unique(unlist(lapply(d, function(x) { unlist(lapply(strsplit(names(x), "_", fixed=TRUE), function(x) {x[1]})) })))
m = matrix(NA, nrow=length(fn), ncol=length(d))
rownames(m) = fn
colnames(m) = sample_id
pb <- txtProgressBar(min=0, max=length(d), style=3)
for (i in 1:length(d)) {
	setTxtProgressBar(pb, i)
	fn = unlist(lapply(strsplit(names(d[[i]]), split="_", fixed=TRUE), function(x) {x[1]}))
 	fx = d[[i]]
 	index = which(duplicated(fn))
 	if (length(index)!=0) {
 		for (j in 1:length(index)) {
 			indx = which(fn==fn[index[j]])
 			fx[indx[which.max(abs(fx[indx]))]] = NA
 		}
 	}
 	fn = fn[!is.na(fx)]
 	fx = fx[!is.na(fx)]
 	m[fn,i] = fx
}
close(pb) 
index = order(apply(m, 1, function(x) {sum(is.na(x))}), decreasing=FALSE)
m = m[index,,drop=FALSE]

index = order(apply(m, 2, function(x) {sum(is.na(x))}), decreasing=FALSE)
m = m[,index,drop=FALSE]
 
for (i in nrow(m):1) {
 	index = order(ifelse(is.na(m[i,]), 1, 0), decreasing=FALSE)
 	m = m[,index,drop=FALSE]
}
 
gene_symbols_env = new.env()
gene_symbols_env$gene_symbols = rownames(m)
 
## partially collapsed mutations
fn = unique(unlist(lapply(d, names)))
mz = matrix(NA, nrow=length(fn), ncol=length(d))
rownames(mz) = fn
colnames(mz) = sample_id
for (i in 1:length(d)) {
 	fx = d[[i]]
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
index = order(apply(mz, 1, function(x) {sum(is.na(x))}), decreasing=FALSE)
mz = mz[index,,drop=FALSE]
 
index = order(apply(mz, 2, function(x) {sum(is.na(x))}), decreasing=FALSE)
mz = mz[,index,drop=FALSE]
 
for (i in nrow(mz):1) {
 	index = order(ifelse(is.na(mz[i,]), 1, 0), decreasing=FALSE)
 	mz = mz[,index,drop=FALSE]
}
 
gene_symbols = gsub(" (3)", "", gsub(" (2)", "", gsub(" (1)", "", rownames(mz), fixed=TRUE), fixed=TRUE), fixed=TRUE)
mt = NULL
for (i in 1:length(gene_symbols_env$gene_symbols)) {
	index = gene_symbols == gene_symbols_env$gene_symbols[i]
	mt = rbind(mt, mz[index,,drop=FALSE])
}

gene_symbols = c("PIK3CA ", "ESR1 ", "PTEN ", "NF1 ", "ARID1A ", "AKT1", "ERBB2", "EGFR")
mz = NULL
for (i in 1:length(gene_symbols)) {
 	index = grep(gene_symbols[i], rownames(mt), fixed=TRUE)
 	mz = rbind(mz, mt[index,,drop=FALSE])
}
mz = rbind(mz, mt[!(rownames(mt) %in% rownames(mz)),,drop=FALSE])
 
for (i in nrow(mz):1) {
 	index = order(ifelse(is.na(mz[i,]), 1, 0), decreasing=FALSE)
 	mz = mz[,index,drop=FALSE]
}

load(url_byl_clinical)
 
tmp = bylc %>%
 	  filter(Nov.ID != "") %>%
 	  select(Nov.ID, Best.Response, Off.Study.Reason, weeks.on.study) %>%
 	  arrange(weeks.on.study)
 	  
mz = mz[,tmp$Nov.ID[tmp$Nov.ID%in%colnames(mz)],drop=FALSE]

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples.pdf"), width=5, height=8)
corr_plot(corr=mz, method = "circle", type = "full", add = FALSE,
        	  col = col_fun(2), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
        	  outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
              addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
              na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=0.90)
dev.off()

col_fun <- colorRampPalette(
	c("#2166AC", "grey95", "#B2182B")
)

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_w_shades.pdf"), width=5, height=8)
corr_plot(corr=mz, border.col = NULL, method = "circle", type = "full", add = FALSE,
          col = "white", bg = NA, title = "", is.corr = FALSE, diag = TRUE,
          outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
          addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
          na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=0.9)
          
export_x = data.frame(mz)
export_x = cbind(gene_symbol = rownames(mz), export_x)
colnames(export_x) = gsub(pattern="X", replacement="", x=colnames(export_x))
write_tsv(export_x, path=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_w_shades.tsv"))
          
mt = mz
colnames(mt) = rep("", ncol(mt))
rownames(mt) = rep("", nrow(mt))
corr_plot(corr=mt, method = "square", type = "full", add = TRUE,
       	  col = col_fun(50), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
       	  outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
          addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
          na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=0.9)
dev.off()

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_w_shades_Key.pdf"), width=5, height=8)
mz["ARAF",1:9] = seq(-4, 4)
corr_plot(corr=mz, method = "square", type = "full", add = FALSE,
       	  col = col_fun(50), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
       	  outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
          addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
          na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=0.9)
dev.off()
