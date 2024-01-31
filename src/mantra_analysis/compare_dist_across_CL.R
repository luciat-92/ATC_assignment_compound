setwd("../../")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(corrplot))
options(stringsAsFactors = F)
options(header = T)

mantra_dist <- fread("data/all_celllines_all_distances.csv", h=T, stringsAsFactors = F, data.table = F)
mantra_dist <- mantra_dist %>% mutate(pair.name = paste(comp1.name, comp2.name, sep = '_and_'))

# count number of pairs per cell line
count_cl <- mantra_dist %>% count(cl.name) %>% arrange(desc(n)) %>%
    mutate(cl.name = factor(cl.name, levels = cl.name))
cl_names <- count_cl$cl.name

pl <- ggplot(count_cl, aes(x = cl.name, y = n))+
    geom_bar(stat = 'identity')+
    scale_y_continuous(trans='log10')+
    theme_classic()+
    theme(axis.title.y = element_blank()) + 
    ylab("N. pairwise distances")+
    coord_flip()
ggsave(pl, file = "data/plot/count_cell_line.png", width = 5, height = 5)

# count number of pairs shared
count_pairs <- mantra_dist %>% count(pair.name)
pl <- ggplot(count_pairs, aes(x = n))+
    geom_bar()+
    scale_y_continuous(trans='log10')+
    theme_classic()+
    ylab("N. unique pairwise distances")+  
    xlab("N. cell lines")
ggsave(pl, file = "data/plot/count_cell_line.png", width = 3, height = 5)


n_cl <- length(cl_names)
common_pair_CL <- matrix(0, nrow = n_cl, ncol = n_cl)
cor_CL <- cor_test_CL <- matrix(0, nrow = n_cl, ncol = n_cl)

for(i in 1:(n_cl-2)){
     print(i)
     for(j in (i+1):n_cl){
         print(j)
          tmp <- mantra_dist %>% filter(cl.name %in% cl_names[c(i, j)]) 
          tmp_int <- tmp %>%
                     filter(pair.name %in% pair.name[duplicated(pair.name)]) %>%
                     group_by(cl.name) %>% arrange(cl.name, pair.name) %>% group_split()
 
          if(length(tmp_int) > 0){
             if(nrow(tmp_int[[1]]) > 1){
                 cor_tmp <- cor.test(tmp_int[[1]][,'d.norm_distance', drop = T], 
                                     tmp_int[[2]][,'d.norm_distance', drop = T],
                                     method = 'spearman')
                 cor_CL[i,j] <- cor_tmp$estimate
                 cor_test_CL[i,j] <- cor_tmp$p.value
            }
            common_pair_CL[i,j] <- nrow(tmp_int[[1]])
        }
    }
}

common_pair_CL <- common_pair_CL + t(common_pair_CL)
colnames(common_pair_CL) <- rownames(common_pair_CL) <- cl_names

cor_CL <- cor_CL + t(cor_CL)
colnames(cor_CL) <- rownames(cor_CL) <- cl_names
colnames(cor_test_CL) <- rownames(cor_test_CL) <- cl_names

png('data/plot/corrplot_distance_spearmancorr.png', width = 1800, height = 1800, res = 200)
corrplot(cor_CL, method = 'square', type = 'upper', order = 'hclust', diag = FALSE, 
        p.mat = cor_test_CL, sig.level = 0.05, pch.cex = 0.8, 
        tl.col = 'black', tl.srt = 45, tl.cex = 0.8)
dev.off()

# save 
obj_corr <- list(cor = cor_CL, pvalue = cor_test_CL, n = common_pair_CL)
save(obj_corr, file = 'data/corr_common_pairs_across_CL.RData')

plot_distance_pairs <- function(cl_name_1, cl_name_2, save_pl = T){

    mantra_dist_subset <- mantra_dist %>% 
                       filter(cl.name %in% c(cl_name_1, cl_name_2)) %>%
                       mutate(cl.name = factor(cl.name, levels = c(cl_name_1, cl_name_2))) %>%
                       filter(pair.name %in% pair.name[duplicated(pair.name)]) %>%
                       group_by(cl.name) %>% arrange(cl.name, pair.name) %>% group_split()
 
    df_subset <- data.frame(cl1 = mantra_dist_subset[[1]][,'d.norm_distance', drop = T], 
                         cl2 = mantra_dist_subset[[2]][,'d.norm_distance', drop = T],
                         pair = mantra_dist_subset[[2]][,'pair.name', drop = T])
 
    pl_dist <- ggplot(df_subset, aes(x = cl1, y = cl2))+
         xlab(sprintf("CL %s", cl_name_1))+ ylab(sprintf("CL %s", cl_name_2))+
         geom_point(size = 0.01) + 
         stat_density2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5)+
         scale_fill_continuous(type = "viridis") +
         ggtitle("Normalized distance")+
         theme_classic() 
 
    if(save_pl){
         ggsave(plot = pl_dist, 
                file = sprintf("data/plot/dist_common_CL%s_CL%s.png", cl_name_1, cl_name_2), 
                height = 5, width = 5)
    }
     
    return(df_subset)
 
}

df <- plot_distance_pairs(cl_name_1 = "independent", cl_name_2 = "cmap-merged")
df <- plot_distance_pairs(cl_name_1 = "independent", cl_name_2 = "A375")
df <- plot_distance_pairs(cl_name_1 = "independent", cl_name_2 = "LNCAP")
