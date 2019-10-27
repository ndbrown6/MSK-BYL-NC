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
			 mutate(oncokb_effect = ifelse(is.na(oncokb_effect), 0, 1)) %>%
			 select(oncokb_effect)
gmaf = cbind(gmaf, oncokb_maf)

pre_ids = unique(cfdna_key$Nov.ID[which(cfdna_key$sample_tx_status == 'Pre' & !cfdna_key$SampleId %in% c(duplicate_id$sample_id, nomut_nopik_x$sample_id))])
post_ids = unique(cfdna_key$Nov.ID[which(cfdna_key$sample_tx_status == 'Post' & !cfdna_key$SampleId %in% c(duplicate_id$sample_id, nomut_nopik_x$sample_id))])

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

## get oncokb `oncogenic' mutations only
for (i in 1:length(d)) {
	index = rep(NA, length(d[[i]]))
	for (j in 1:length(d[[i]])) {
		index[j] = gmaf$oncokb_effect[which(paste0(gmaf$Hugo_Symbo, "_", gmaf$HGVSp_Short) %in% names(d[[i]][j]))[1]]
	}
	d[[i]] = d[[i]][index==1]
}
 
## fully collapsed mutations
fn = unique(unlist(lapply(d, function(x) { unlist(lapply(strsplit(names(x), "_", fixed=TRUE), function(x) {x[1]})) })))
m = matrix(NA, nrow=length(fn), ncol=length(d))
rownames(m) = fn
colnames(m) = sample_id
for (i in 1:length(d)) {
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
uuid_2 = c("PIK3CA", "PTEN", "ARID1A")
uuid_3 = c("ESR1", "NF1")
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
rownames(tmp.1) = paste0(rep(uuid_2, each=2), " (", rep(1:2, times=length(uuid_2)), ")")
 
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
rownames(tmp.2) = paste0(rep(uuid_3, each=3), " (", rep(1:3, times=length(uuid_3)), ")")
 
mz = rbind(tmp.1, tmp.2, tmp.0)
colnames(mz) = sample_id
index = order(apply(mz, 1, function(x) {sum(is.na(x))}), decreasing=FALSE)
mz = mz[index,,drop=FALSE]
 
for (i in nrow(mz):1) {
 	index = order(ifelse(is.na(mz[i,]), 1, 0), decreasing=FALSE)
 	mz = mz[,index,drop=FALSE]
}

mt = mz
 
gene_symbols = c("PIK3CA ",
				 "ESR1 ",
				 "PTEN ",
				 "NF1 ",
				 "ARID1A ",
				 "AKT1",
				 "ERBB2",
				 "EGFR")

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

index = order(apply(mz, 1, function(x) { sum(x<0, na.rm=TRUE) }), decreasing=FALSE)
mz = mz[index,,drop=FALSE]
index = order(apply(mz, 1, function(x) { sum(x>=0, na.rm=TRUE) }), decreasing=TRUE)
mz = mz[index,,drop=FALSE]
mz = mz[c(rownames(mz)[rownames(mz)!="TP53"],"TP53"),,drop=FALSE]

load(url_byl_clinical)
 
tmp = bylc %>%
 	  filter(Nov.ID != "") %>%
 	  select(Nov.ID, Best.Response, Off.Study.Reason, weeks.on.study) %>%
 	  arrange(weeks.on.study)
 	  
mz = mz[,tmp$Nov.ID[tmp$Nov.ID%in%colnames(mz)],drop=FALSE]
index = rownames(mz) %in% c("CDH1", "BRCA1", "BRCA2", "ATM")
mz = mz[!index,,drop=FALSE]

gene_symbols = c("PIK3CA (1)",
				 "PIK3CA (2)",
				 "ESR1 (1)",
				 "ESR1 (2)",
				 "ESR1 (3)",
				 "PTEN (1)",
				 "PTEN (2)",
				 "AKT1",
				 "MTOR",
				 "BRAF",
				 "HRAS",
				 "MAP2K1",
				 "NF1 (1)",
				 "NF1 (2)",
				 "NF1 (3)",
				 "FGFR1",
				 "ARID1A (1)",
				 "ARID1A (2)")

mz = rbind(mz[gene_symbols,,drop=FALSE], mz[!(rownames(mz) %in% gene_symbols),,drop=FALSE])

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_oncokb.pdf"), width=5, height=5.5)
corr_plot(corr=mz, method = "circle", type = "full", add = FALSE,
        	  col = col_fun(2), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
        	  outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
              addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
              na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=1)
dev.off()

col_fun <- colorRampPalette(
	c("#3288bd", "white", "#d53e4f")
)

sample_names = c("045", "006", "043",
				 "060", "026", "031",
				 "032", "033", "047",
				 "059", "008", "053",
				 "054", "002", "057",
				 "050", "040", "016",
				 "023", "025", "021",
				 "035", "038", "007",
				 "034", "009", "046",
				 "039", "022", "029",
				 "048", "004")
mz = mz[,sample_names,drop=FALSE]

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_w_shades_oncokb.pdf"), width=5, height=5.5)
corr_plot(corr=mz, border.col = NULL, method = "circle", type = "full", add = FALSE,
          col = "white", bg = NA, title = "", is.corr = FALSE, diag = TRUE,
          outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
          addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
          na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=1)
          
export_x = data.frame(mz)
export_x = cbind(gene_symbol = rownames(mz), export_x)
colnames(export_x) = gsub(pattern="X", replacement="", x=colnames(export_x))
write_tsv(export_x, path=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_w_shades_oncokb.tsv"))          
          
mt = mz
colnames(mt) = rep("", ncol(mt))
rownames(mt) = rep("", nrow(mt))
corr_plot(corr=mt, method = "square", type = "full", add = TRUE,
       	  col = col_fun(50), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
       	  outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
          addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
          na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=1)
dev.off()

pdf(file=str_c(out_dir, "VAF_matrix_by_sample_uncollapsed_all_samples_w_shades_oncokb_Key.pdf"), width=5, height=5.5)
mz["TP53",1:9] = seq(-4, 4)
corr_plot(corr=mz, method = "square", type = "full", add = FALSE,
       	  col = col_fun(50), bg = NA, title = "", is.corr = FALSE, diag = TRUE,
       	  outline = TRUE, mar = c(0, 0, 0, 0), addgrid.col = "grey90",
          addCoef.col = NULL, addCoefasPercent = FALSE, order = "original",
          na.label = " ", cl.pos="n", tl.cex=.75, tl.col="black", win.asp=1)
dev.off()
