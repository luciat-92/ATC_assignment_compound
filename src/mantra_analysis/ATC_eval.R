#setwd("../../")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(umap))
options(stringsAsFactors = F)
options(header = T)

### functions ###
get_atc <- function(compound, df){

    lev5 <- df %>% filter(mantra_name %in% compound) %>% 
    filter(!duplicated(level5)) %>% select(level5) 
    atc_codes <- data.frame(compound = compound, 
                            level5 = paste0(lev5$level5[lev5$level5 != ""], collapse = "-"), 
                            n_atc=sum(lev5$level5 != ""))
    return(atc_codes)

}

check_same_atc <- function(code1, code2, depth){
 
    if(nchar(code1) > 7){
        code1 <- strsplit(code1, split = '-')[[1]]
    }
    if(nchar(code2) > 7){
        code2 <- strsplit(code2, split = '-')[[1]]
    }

    if(length(code1) == 1 & length(code2) == 1){
        condition <- substring(code1,1,depth) == substring(code2,1,depth)
    }else if (length(code1) > 1 & length(code2) == 1){
       condition <- sapply(code1, function(x) substring(x,1,depth) == substring(code2,1,depth))
    }else if (length(code1) == 1 & length(code2) > 1){
        condition <- sapply(code2, function(x) substring(code1,1,depth) == substring(x,1,depth))
    }else{
        condition <- sapply(code1, function(x) sapply(code2, function(y) substring(x,1,depth) == substring(y,1,depth)))
    }

    return(condition)

}

compute_PPV <- function(atc_match_sorted, n_closest){

    df <- data.frame(PPV = sapply(1:n_closest, function(n) sum(atc_match_sorted[1:n])/n), 
                     k = 1:n_closest)

    return(df)
}

get_compounds_with_atcs <- function(mantra_dist, df_chembl, cl_name = NULL){
    
    if(!is.null(cl_name)){
        mantra_dist <- mantra_dist %>% filter(cl.name == cl_name)
    }

    compound_unique <- c(mantra_dist$comp1.name, mantra_dist$comp2.name) %>% unique()
    df_chembl <- df_chembl %>% filter(mantra_name != "")
    common_c <- intersect(df_chembl$mantra_name, compound_unique)

    if(length(common_c) == 0){
        mantra_atc <- data.frame(compound = c(), level5 = c(), n_atc = c(), cl_name = c())
    }else{
        mantra_atc <- lapply(common_c, function(x) get_atc(x, df_chembl))
        mantra_atc <- do.call(rbind, mantra_atc)
        if(!is.null(cl_name)){
            mantra_atc$cl_name <- cl_name
       }else{
            mantra_atc$cl_name <- "All cell lines"
        }
    }
    
    return(mantra_atc)

}

create_dist_matrix <- function(mantra_df, cl_name){

    mantra_df <- mantra_df %>% filter(cl.name == cl_name)
    compound_names <- unique(c(mantra_df$comp1.name, mantra_df$comp2.name))
    n_comp <- length(compound_names)
    dist_matrix <- diag(x = 0, nrow = n_comp)
    colnames(dist_matrix) <- rownames(dist_matrix) <- compound_names
    for(name_c in compound_names){
        comp1_tmp <- mantra_df[mantra_df$comp1.name == name_c,]
        if(nrow(comp1_tmp)>0){
        dist_matrix[name_c, match(comp1_tmp$comp2.name, colnames(dist_matrix))] <- comp1_tmp$d.norm_distance
        }
    }

    dist_matrix <- dist_matrix + t(dist_matrix)
    if(!isSymmetric(dist_matrix)){stop("distance matrix not symmetric")}
    return(dist_matrix)

}

#######################

df_chembl <- fread("data/chembl_atc_mapping.tsv", h=T, stringsAsFactors = F, data.table = F)
mantra_dist <- fread("data/all_celllines_all_distances.csv", h=T, stringsAsFactors = F, data.table = F)

cl_names <- unique(mantra_dist$cl.name)

mantra_atc <- get_compounds_with_atcs(mantra_dist, df_chembl)
mantra_atc_CL <- lapply(cl_names, function(x) get_compounds_with_atcs(mantra_dist, df_chembl, cl_name = x))
mantra_atc_CL <- do.call(rbind, mantra_atc_CL)
mantra_atc_CL <- mantra_atc_CL %>% filter(n_atc >0)
df_count <- mantra_atc_CL %>% count(cl_name) %>% 
                arrange(desc(n))
df_count$cl_name <- factor(df_count$cl_name, levels = df_count$cl_name)

pl <- ggplot(df_count, aes(x = cl_name, y = n))+
    geom_bar(stat = "identity")+
    theme_classic()+
    theme(axis.title.y = element_blank()) + 
    ylab("N. compounds with ATC")+
    coord_flip()
ggsave(pl, file = "data/plot/count_compounds_with_atc.png", width = 5, height = 5)

# add ATC info on cell line, remove pairs without ATC info
mantra_dist_withATC <- mantra_dist %>% 
    mutate(comp1.atc = mantra_atc$level5[match(mantra_dist$comp1.name,mantra_atc$compound)]) %>%
    mutate(comp2.atc = mantra_atc$level5[match(mantra_dist$comp2.name,mantra_atc$compound)]) %>%
    mutate(comp1.atc = case_when(is.na(comp1.atc) ~ "", TRUE ~ as.character(comp1.atc))) %>%
    mutate(comp2.atc = case_when(is.na(comp2.atc) ~ "", TRUE ~ as.character(comp2.atc))) %>%
    filter(comp1.atc != ""  & comp2.atc != "" )

# add T if atc codes match for a certain length
mantra_dist_withATC$atc_1_length <- sapply(1:nrow(mantra_dist_withATC), function(x) 
                                            any(check_same_atc(mantra_dist_withATC$comp1.atc[x], 
                                                               mantra_dist_withATC$comp2.atc[x], 1)))

mantra_dist_withATC$atc_3_length <- sapply(1:nrow(mantra_dist_withATC), function(x) 
                                            any(check_same_atc(mantra_dist_withATC$comp1.atc[x], 
                                                               mantra_dist_withATC$comp2.atc[x], 3)))

mantra_dist_withATC$atc_4_length <- sapply(1:nrow(mantra_dist_withATC), function(x) 
                                            any(check_same_atc(mantra_dist_withATC$comp1.atc[x], 
                                                               mantra_dist_withATC$comp2.atc[x], 4)))

mantra_dist_withATC$atc_5_length <- sapply(1:nrow(mantra_dist_withATC), function(x) 
                                            any(check_same_atc(mantra_dist_withATC$comp1.atc[x], 
                                                               mantra_dist_withATC$comp2.atc[x], 5)))                                                               
                                                               
# sort based on distance
mantra_dist_withATC <- mantra_dist_withATC %>% arrange(cl.name, d.norm_distance)

# plot PPV for each cell line
n_closest <- 20000
depth_atc <- c('atc_1_length','atc_3_length', 'atc_4_length','atc_5_length')
cl_names <- unique(mantra_dist_withATC$cl.name)

df_PPV <- data.frame(PPV = c(), k = c(), type = c(), cl = c())
for(id_cl in 1:length(cl_names)){
    print(cl_names[id_cl])
    for(col in depth_atc){
        tmp <- compute_PPV(mantra_dist_withATC[mantra_dist_withATC$cl.name == cl_names[id_cl],col], n_closest)
        tmp$type <- col
        tmp$cl <- cl_names[id_cl]
        df_PPV <- rbind(df_PPV, tmp)
    }
}

df_PPV$type <- factor(df_PPV$type, levels = depth_atc)
df_PPV$cl <- factor(df_PPV$cl, levels = cl_names)

# plot
pl <- ggplot(df_PPV, aes(x = k, y = PPV, color = type)) +
    geom_line()+
    facet_wrap(.~cl, nrow = 4, ncol=8)+
    theme_bw() + theme(legend.position = 'bottom') + 
    xlab('n. of connections') +
    ylab('Positive Predicted Values')
ggsave(plot = pl, file = 'data/plot/PPV_total_atc_match_allCL.png', height = 10, width = 15)

pl_zoom <- filter(df_PPV, k <=100) %>%  
    ggplot(aes(x = k, y = PPV, color = type)) +
    geom_line()+
    facet_wrap(.~cl, nrow = 4, ncol=8)+
    theme_bw() + theme(legend.position = 'bottom') + 
    xlab('n. of connections') +
    ylab('Positive Predicted Values')
ggsave(plot = pl_zoom, file = 'data/plot/PPV_zoom100_atc_match_allCL.png', height = 10, width = 15)


#### find pairs that are in cmap-merged but not cell line idependent 
mantra_dist_withATC <- mantra_dist_withATC %>% mutate(pair.name =paste(comp1.name, comp2.name, sep = "_and_"))

mantra_dist_withATC_cmap <- mantra_dist_withATC %>% filter(cl.name == 'cmap-merged')
mantra_dist_withATC_ind <- mantra_dist_withATC %>% filter(cl.name == 'independent')
diff_merge <- anti_join(mantra_dist_withATC_cmap, mantra_dist_withATC_ind, by = "pair.name")
common_merge <- inner_join(mantra_dist_withATC_cmap, mantra_dist_withATC_ind, by = "pair.name")

# add diff merge to mantra_dist_withATC_ind
new_merged <- rbind(mantra_dist_withATC_ind, diff_merge) %>% arrange(d.norm_distance)
df_PPV_merged <- data.frame(PPV = c(), k = c(), type = c())
for(col in depth_atc){
        tmp <- compute_PPV(new_merged[,col], n_closest)
        tmp$type <- col
        df_PPV_merged <- rbind(df_PPV_merged, tmp)
}
df_PPV_merged$type <- factor(df_PPV_merged$type, levels = depth_atc)

pl_zoom <- filter(df_PPV_merged, k <=400) %>%  
    ggplot(aes(x = k, y = PPV, color = type)) +
    geom_line()+
    theme_bw() + theme(legend.position = 'bottom') + 
    xlab('n. of connections') +
    ylab('Positive Predicted Values')
ggsave(plot = pl_zoom, file = 'data/plot/PPV_zoom100_atc_match_merged.png', height = 5, width = 6)

# compare top cmap in independent
top_cmap <- mantra_dist_withATC_cmap[1:14,]
top_cmap_in_ind <-  mantra_dist_withATC_ind %>% filter(pair.name %in% top_cmap$pair.name)
write.table(top_cmap, file = 'data/top_cmap_withATC.txt',col.names = T, row.names = F, sep = '\t', quote = F)

top_ind <- mantra_dist_withATC_ind[1:14,]
top_ind_in_cmap <-  mantra_dist_withATC_cmap %>% filter(pair.name %in% top_ind$pair.name)
write.table(top_ind, file = 'data/top_ind_withATC.txt',col.names = T, row.names = F, sep = '\t', quote = F)

# plot PPV for common paris
common_pairs <- intersect(mantra_dist_withATC_cmap$pair.name, mantra_dist_withATC_ind$pair.name)
mantra_dist_withATC_cmap_f <- mantra_dist_withATC_cmap %>% filter(pair.name %in% common_pairs)
mantra_dist_withATC_ind_f <- mantra_dist_withATC_ind %>% filter(pair.name %in% common_pairs)

df_PPV_cmap_f <- data.frame(PPV = c(), k = c(), type = c())
for(col in depth_atc){
        tmp <- compute_PPV(mantra_dist_withATC_cmap_f[,col], n_closest)
        tmp$type <- col
        df_PPV_cmap_f <- rbind(df_PPV_cmap_f, tmp)
}
df_PPV_cmap_f$type <- factor(df_PPV_cmap_f$type, levels = depth_atc)
df_PPV_cmap_f$cl_name <- "cmap-merged"

df_PPV_ind_f <- data.frame(PPV = c(), k = c(), type = c())
for(col in depth_atc){
        tmp <- compute_PPV(mantra_dist_withATC_ind_f[,col], n_closest)
        tmp$type <- col
        df_PPV_ind_f <- rbind(df_PPV_ind_f, tmp)
}
df_PPV_ind_f$type <- factor(df_PPV_ind_f$type, levels = depth_atc)
df_PPV_ind_f$cl_name <- "independent"
df_PPV_f <- rbind(df_PPV_ind_f, df_PPV_cmap_f)

pl <- ggplot(df_PPV_f, aes(x = k, y = PPV, color = type)) +
    geom_line()+
    facet_wrap(.~cl_name, nrow = 1)+
    theme_bw() + theme(legend.position = 'bottom') + 
    xlab('n. of connections') +
    ylab('Positive Predicted Values')
ggsave(plot = pl, file = 'data/plot/PPV_atc_match_common.png', height = 5, width = 6)

pl_zoom <- filter(df_PPV_f, k <=100) %>%  
    ggplot(aes(x = k, y = PPV, color = type)) +
    geom_line()+
    facet_wrap(.~cl_name, nrow = 1)+
    theme_bw() + theme(legend.position = 'bottom') + 
    xlab('n. of connections') +
    ylab('Positive Predicted Values')
ggsave(plot = pl_zoom, file = 'data/plot/PPV_zoom100_atc_match_common.png', height = 5, width = 6)

#################
### UMAP plot ###
#################

custom.config <- umap.defaults
custom.config$n_neighbors <- 10
custom.config$min_dist <- 0.05

dist_cmap_merged <- create_dist_matrix(mantra_df = mantra_dist_withATC, cl_name = 'cmap-merged')
dist_independent <- create_dist_matrix(mantra_df = mantra_dist_withATC, cl_name = 'independent')

atc_cmap_merged <- mantra_atc_CL %>% filter(cl_name == 'cmap-merged')
atc_cmap_merged <- atc_cmap_merged[match(colnames(dist_cmap_merged),atc_cmap_merged$compound), ]

atc_independent <- mantra_atc_CL %>% filter(cl_name == 'independent')
atc_independent <- atc_independent[match(colnames(dist_independent),atc_independent$compound), ]

umap_cmap_merged <- umap(dist_cmap_merged, config=custom.config, input="dist")
umap_independent <- umap(dist_independent, config=custom.config, input="dist")

df_umap_cmap_merged <- data.frame(UMAP_1 = umap_cmap_merged$layout[,1],
                      UMAP_2 = umap_cmap_merged$layout[,2],
                      atc_1_length = substring(atc_cmap_merged$level5, 1,1))
df_umap_cmap_merged$atc_1_length <- factor(df_umap_cmap_merged$atc_1_length)

df_umap_independent <- data.frame(UMAP_1 = umap_independent$layout[,1],
                      UMAP_2 = umap_independent$layout[,2],
                      atc_1_length = substring(atc_independent$level5, 1,1))
df_umap_independent$atc_1_length <- factor(df_umap_independent$atc_1_length)

pl_cmap <- ggplot(df_umap_cmap_merged, aes(x = UMAP_1, y = UMAP_2, color = atc_1_length)) +
    geom_point(size=0.2)+
    theme_bw() +
    xlab('UMAP 1') + ylab('UMAP 2') +  ggtitle("CL cmap-merged")

pl_ind <- ggplot(df_umap_independent, aes(x = UMAP_1, y = UMAP_2, color = atc_1_length)) +
    geom_point(size=0.2)+
    theme_bw() +
    xlab('UMAP 1') + ylab('UMAP 2') + ggtitle("CL independent")

pl <- ggarrange(plotlist = list(pl_ind, pl_cmap), common.legend = T, nrow = 1)
ggsave(plot = pl, file = 'data/plot/UMAP_norm_dist_ATC1_colored.png', height = 6, width = 10)
