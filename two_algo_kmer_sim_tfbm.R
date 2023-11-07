# K-mer similarity------------
# two algorithms in your imp K-mer list
library(dplyr)
library(tibble)
library(universalmotif)
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridtext)


dir <- '/path/to/all/folder/of/imp/file/group/'
file_name <- 'group_folder_name'
filenames <- list.files(paste0(dir,file_name),
                        pattern="*.txt",
                        full.names=FALSE)
head(filenames)

# rank K-mer importance in each different algorithms
df1 <- read.delim(paste0(dir, filenames[1])) %>%
  select(1,4)
df2 <- read.delim(paste0(dir, filenames[2])) %>%
  select(1,4)
names(df1) <- c('motif', 'median1')
names(df2) <- c('motif', 'median2')


df <- merge(df1,df2, by = 'motif')

df$rank1 <- order(df$median1)
df$rank2 <- order(df$median2)
df$rank_mean <- apply(df[,c('rank1', 'rank2')], 1, mean)
df <- df[order(df$rank_mean),]
df <- df$motif
# option: df <- df$motif[1:10,] as top10 K-mers
# motif pwm can be apply to generate consensus motif in top10 K-mers


# TFBM database from ShuiLab (Arabidopsis)
tfbm <- read_meme('/full/path/to/Athaliana_TFBM_v1.01.tm.index.direct.index_modified_meme.txt')

# function compare_motifs() included multiple parameters.
# please check the manuals for comprehensive description.  
comparisons <- compare_motifs(c(list.motif, tfbm),
                              method = "PCC", # Pearson Correlation Coefficient
                              min.mean.ic = 0,
                              score.strat = "a.mean")


# remove K-mer to K-mer PCC
comparisons2 <- comparisons[-c(1:length(df)), ]


for (i in 1:length(df)) {
  if(i == 1){
    
    df <- data.frame(comparisons2[,i])
    cname <- list.motif[[i]]@name
    colnames(df) <- cname
    rname <- row.names(df)
    df$motif <- rname
    df <- df[order(-df[cname]), ][1,] # top 10 similar motifs
    df <- df[,c('motif', cname)]
    row.names(df) <- 1
    
  } else {
    
    df2 <- data.frame(comparisons2[,i])
    cname2 <- list.motif[[i]]@name
    colnames(df2) <- cname2
    rname2 <- row.names(df2)
    df2$motif <- rname2
    df2 <- df2[order(-df2[cname2]), ][1,]
    df2 <- df2[,c('motif', cname2)]
    row.names(df2) <- 1
    
    df <- full_join(df, df2, by = 'motif')
  }
  
}


df3 <- gather(df, key = 'kmers', value = 'PCC', -motif)
df3 <- df3 %>% filter(PCC != '')
df3$cluster <- file_name


# TFBM motif family----
motif.family <- df3

string <- motif.family$motif
# remove unwanted dash-line and any digit behind it as family name
string <- str_replace_all(string, '\\_m1', '')
head(string)

motif.family$motif <- string

motif.family <- motif.family %>%
  group_by(kmers, motif) %>%
  distinct() %>% 
  summarise_all(median)

motif.family <- motif.family %>% 
  column_to_rownames(var = 'kmers')

output <- data.frame(motif.family)
output <- rownames_to_column(output, var = 'k-mers')

write.table(output,
            paste0(dir, file_name, '_top_sim_two_algo_tfbm.txt'),
            row.names = F,
            quote = F,
            sep = '\t')

# tfbm motif family heatmap-------

col_fun = colorRamp2(c(0, 0.5, 1), c('blue', 'white', 'red'))
ha = rowAnnotation(PCC = as.numeric(motif.family$PCC),
                   col = list(PCC = col_fun),
                   annotation_legend_param = list(
                     PCC = list(direction = "horizontal")
                   ))


mfh <- data.frame(motif.family[,1])
colnames(mfh) <- 'motif'
row.names(mfh) <- rownames(motif.family)


ht <- Heatmap(mfh,
              name = 'family',
              column_title = gt_render(
                paste0("<span style = 'font-size:10pt'>**TFBM**</span><br>")),
              show_column_names = FALSE,
              row_title = 'Important Kmers',
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 8),
              show_row_names = TRUE,
              show_heatmap_legend = FALSE,
              row_order = rownames(motif.family),
              right_annotation = ha,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(mfh[i, j]), x, y, gp = gpar(fontsize = 10))
              }
)


pdf(paste0(dir, file_name, '_top_sim_tow_algo_tfbm.pdf'))
#option: add type="cairo" needed in HPC system
draw(ht,
     merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = 'bottom')

dev.off()

save.image(file = paste0(dir, paste(file_name, 'tow_algo_tfbm.RData', sep = '_')))







