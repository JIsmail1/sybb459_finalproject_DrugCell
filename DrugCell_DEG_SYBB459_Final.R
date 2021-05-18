# Makaela Mews
# SYBB459  Final
# 5/18/21

# Download DrugCell Info: https://github.com/idekerlab/DrugCell
# Download Affy pre-processed RMA normalised expression data: https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html
# Download Affy cell line names from E-MTAB-3610.sdrf.txt found at: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3610/ 

# Set up environment
setwd('/home/jovyan/SYBB459')
install.packages('tidyr')
install.packages('dplyr')
library(tidyr)
library(dplyr)
install.packages('sqldf')
library(sqldf)
install.packages('stringr')
library(stringr)

# Load Affy data
cell <- read.delim(file ='Cell_line_RMA_proc_basalExp_ref.txt', check.names = FALSE)

# used check.names = FALSE so that IDs match key
cell$GENE_SYMBOLS <- toupper(cell$GENE_SYMBOLS)

# Identify cell lines that contain duplicates i.e. #.1 
cell_ex_id <- data.frame(str_split_fixed(colnames(cell), "[.]", n = 2))
head(cell_ex_id)
dim(cell_ex_id)

# subset of dataframe where there is duplicates
cell_dup <- cell_ex_id[which(cell_ex_id$X2 == '1'), ]
cell_dup
#rownames = colnames of cell so can subset the columns out
dim(cell)

# Affy data de-duplicated
cell_v2 <- cell[,-c(as.numeric(rownames(cell_dup)))]
#head(cell_v2)
dim(cell_v2)


# Load DrugCell gene index data
gene2ind <- read.delim(file = 'gene2ind.txt', header = FALSE)

# Load Drug Cell cell index + tissue type
cell2ind <- read.delim(file = 'cell2ind.txt', header = FALSE)
colnames(cell2ind) <- c('index' , 'cell_line')

cell_mut <- read.csv(file ='cell2mutation.txt', check.names = FALSE, header = FALSE)

# Formating MUT binary vec
cell_mut <- read.csv(file ='cell2mutation.txt', check.names = FALSE, header = FALSE)

mut_vec_names <- gene2ind$V1

colnames(cell_mut) <- mut_vec_names
rownames(cell_mut) <- cell2ind$index
write.csv(cell_mut, file = 'cell_mut2.csv', quote = FALSE, row.names = TRUE)

# Formatting of MUT binary vec, more human readable
cell_mut_R <- cell_mut
mut_vec_names2 <- gene2ind$V2
colnames(cell_mut_R) <- mut_vec_names2
rownames(cell_mut_R) <- cell2ind$cell_line


#Used Linux to pre-process files 
#awk -F '\t' 'OFS="\t" '{ print $1, $3}' E-MTAB-3610.sdrf.txt > E-MTAB-3610.cell2.txt
#sed 's/-//g' E-MTAB-3610.cell2.txt > E-MTAB-3610.cellref2.txt

#awk -F '\t' '{ print $2 }' E-MTAB-3610.cellref2.txt > E-MTAB-3610.cellname2.txt

#awk -F '\t' '{ print $1 }' E-MTAB-3610.cellref2.txt > E-MTAB-3610.cellID2.txt
#sed 's/ /_/g' E-MTAB-3610.cellID2.txt > E-MTAB-3610.cellID3.txt


#awk '{ print $1 }' E-MTAB-3610.cellID3.txt | sed 's/_/ /g' | awk '{ print $NF}' > E-MTAB-3610.cellID4.txt

#less -N E-MTAB-3610.cellID4.txt

#wc -l E-MTAB-3610.cellID4.txt
#1019 E-MTAB-3610.cellID4.txt

#wc -l E-MTAB-3610.cellname2.txt
#1019 E-MTAB-3610.cellname2.txt

# 872 J132_EPH10P5_SKMEL28_Skin_905954        SKMEL28
# 995 J154_EPH04P12_SKMEL28_skin_905954       SKMEL28

# Load tissue / cell-line files 
# Issue affy data IDs do NOT directly correspond with DrugCell identifers
# only portion of the full E-MTAB name matches between
# i.e cell2ind.txt file 

dcell_id <- read.delim(file = 'E-MTAB-3610.cellID4.txt', header = TRUE)
dcell_name <- read.delim(file = 'E-MTAB-3610.cellname2.txt', header = TRUE)
dcell_name$Characteristics.cell.line. <- toupper(dcell_name$Characteristics.cell.line.)
#glimpse(dcell_id)
#glimpse(dcell_name)

dcell <- cbind(dcell_id, dcell_name)
colnames(dcell) <- c('ID' , 'cell_line')

# Split cell_line from tissue
dcell_line_tissue <- data.frame(str_split_fixed(cell2ind$cell_line, "_", n = 2))

dcell_ind <- cbind(cell2ind$index, dcell_line_tissue)
colnames(dcell_ind) <- c('dindex', 'cell_line', 'tissue')

colnames(dcell) <- c('ID' , 'cell_line')

dcell_ict <- sqldf("SELECT dcell.ID, dcell_ind.dindex, dcell_ind.cell_line, dcell_ind.tissue FROM dcell INNER JOIN dcell_ind ON dcell.cell_line = dcell_ind.cell_line")
dim(dcell_ict)

#removed duplicate from Dcell data so that I can have consistency between affy array  & dcell
dcell_ict_nodup <- dcell_ict[!duplicated(dcell_ict),]

t_cell <- data.frame(t(cell_v2))
t_cell$Sample <- rownames(t_cell)
head(t_cell)
#glimpse(t_cell)

t_cell_ID <- data.frame(t_cell$Sample)
colnames(t_cell_ID) <- c('Sample')

#glimpse(t_cell_ID)

#dcell_ict_nodup
top <- t_cell[1:2, ]
head(top)

# Dataframe containing keys between colname(cell_v2) to correspond with cell2ind.txt (DrugCell) file
merge <- sqldf("SELECT * FROM t_cell_ID INNER JOIN dcell_ict_nodup ON t_cell_ID.Sample = dcell_ict_nodup.ID ")
head(merge)

# retains all entries so that creation on binary vector is more straight forward 
suc <- sqldf("SELECT * FROM gene2ind LEFT JOIN cell_v2 ON gene2ind.V2 = cell_v2.GENE_SYMBOLS")




# Access tissue frequency distribution, preliminary for baseline creation 
ttab <- data.frame(table(dcell_ict_nodup$tissue))
dim(ttab)

#Tissue with sample size at least 5
ttab2 <- ttab[which(ttab$Freq >= 5),]
ttab2

ttab3 <- ttab[which(ttab$Freq >= 10),]
ttab3

# If we restrict the analysis to include only tissues that contain at least 5 different cell lines 
sum(ttab2$Freq)
# 939 new sample size

# If we restrict the analysis to include only tissues that contain at least 10 different cell lines 
sum((ttab[which(ttab$Freq >= 10),])$Freq)
#916 new sample size

#affy data with only sample information not gene info 
suc_sample <- suc[,-c(1:4)]

# gene info for affy data
suc_geneinfo <- suc[,c(1:2)]
head(suc_sample)

colnames(suc_geneinfo) <- c('Dcell_gene_indx', 'Dcell_hgnc_symbol')
dim(suc_sample)



#Convert colnames from affy expression data to drugcell data
suc_samt <- data.frame(t(suc_sample[,]),check.names = FALSE)
aff_id <- data.frame(as.integer(rownames(suc_samt)))
colnames(aff_id) <- 'aff_id'
suc_samt2 <- cbind(aff_id,suc_samt)
head(suc_samt)
head(suc_samt2)
suc_samt3_id <- suc_samt2[,1:2]

suc_sam_fi <- sqldf("SELECT * FROM suc_samt3_id LEFT JOIN merge ON suc_samt3_id.aff_id = merge.ID ")
head(suc_samt3_id)

head(suc_sam_fi)
dim(suc_sam_fi)

colnames(suc_sample) <- suc_sam_fi$dindex
head(suc_sample)

# Vector of tissues that contain at least 10 dif. cell lines
tiss_vec <- as.character(ttab3$Var1)
# Key df
head(merge)

my.sign <- function(x) {
ifelse(x > 2, 1, ifelse(x < -2, 1, ifelse(x, 0)))}


for (i in tiss_vec){
    # Vector of colnames to extract from Affy data for 1 tissue
    vec_col_aff <- merge[which(merge$tissue == i), 3]
    tissue <- i
    s_vec <- NULL
    
    for (i in vec_col_aff){
    string <- toString(i)
    s_vec <- append(s_vec, string)
    }
    s_vec
    suc_tiss <- suc_sample %>% select(s_vec)
    suc_mat <- t(suc_tiss)
    suc_mat <- scale(suc_mat, center = TRUE, scale = TRUE)
    suc_mat2 <- apply(suc_mat, c(1, 2), my.sign )    # Apply own function to each element
    suc_mat2[is.na(suc_mat2)] <- 0
    suc_done <- data.frame(t(suc_mat2), check.names = FALSE)
    #suc_done <- cbind(suc_geneinfo, suc_ref)
    name <- paste0("Tissue3_", tissue)
    write.csv(suc_done, file = name, quote = FALSE, row.names = FALSE)
    
}

df1 <- read.csv(file = 'Tissue3_AUTONOMIC_GANGLIA', check.names = FALSE)
df2 <- read.csv(file = 'Tissue3_BONE', check.names = FALSE)
df3 <- read.csv(file = 'Tissue3_BREAST', check.names = FALSE)
df4 <- read.csv(file = 'Tissue3_CENTRAL_NERVOUS_SYSTEM', check.names = FALSE)
df5 <- read.csv(file = 'Tissue3_CERVIX', check.names = FALSE)
df6 <- read.csv(file = 'Tissue3_ENDOMETRIUM', check.names = FALSE)
df7 <- read.csv(file = 'Tissue3_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', check.names = FALSE)
df8 <- read.csv(file = 'Tissue3_KIDNEY', check.names = FALSE)
df9 <- read.csv(file = 'Tissue3_LARGE_INTESTINE', check.names = FALSE)
df10 <- read.csv(file = 'Tissue3_LIVER', check.names = FALSE)
df11 <- read.csv(file = 'Tissue3_LUNG', check.names = FALSE)
df12 <- read.csv(file = 'Tissue3_OESOPHAGUS', check.names = FALSE)
df13 <- read.csv(file = 'Tissue3_OVARY', check.names = FALSE)
df14 <- read.csv(file = 'Tissue3_PANCREAS', check.names = FALSE)
df15 <- read.csv(file = 'Tissue3_SKIN', check.names = FALSE)
df16 <- read.csv(file = 'Tissue3_SOFT_TISSUE', check.names = FALSE)
df17 <- read.csv(file = 'Tissue3_STOMACH', check.names = FALSE)
df18 <- read.csv(file = 'Tissue3_THYROID', check.names = FALSE)
df19 <- read.csv(file = 'Tissue3_UPPER_AERODIGESTIVE_TRACT', check.names = FALSE)
df20 <- read.csv(file = 'Tissue3_URINARY_TRACT', check.names = FALSE)

df_fi <- cbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df17, df18, df19, df20)
dim(df_fi)
#glimpse(df_fi)

df_fi_nodup <- df_fi[!duplicated(as.list(df_fi))]
head(df_fi_nodup)
dim(df_fi_nodup)

#my.signx <- function(x) {
#ifelse(x > 1, 1, ifelse(x == 1, 1, ifelse(x, 0)))}


#my.sign <- function(x) {
#ifelse(x > 2, 1, ifelse(x < -2, 1, ifelse(x, 0)))}


mut_vec_names <- gene2ind$V1
colnames(cell_mut) <- mut_vec_names
dfs <- data.frame(t(mut_vec_names), check.names = FALSE)
dfs_af <- dfs
dfs_dc <- dfs
col <- colnames(df_fi_nodup)
cell_mut2 <- data.frame(t(cell_mut), check.names = FALSE)
#    suc_ref <- data.frame(t(suc_mat2), check.names = FALSE)
#    suc_done <- cbind(suc_geneinfo, suc_ref)
#data1[data1 == "A"] <- "XXX"

for (i in col){
    sel_col_af <- df_fi_nodup %>% select(i)
    #cell line (col) x genes (row)
    sel_col_dc <- cell_mut2 %>% select(i)
    sel_combo <- cbind(sel_col_dc, sel_col_af)
    sel_combo$exp <- apply(sel_combo, 1, sum)
    sel_combo[sel_combo == 2] <- 1
    sel_mat3 <- sel_combo[,3]
    sel_mat4 <- data.frame(t(sel_mat3), check.names = FALSE)
    rownames(sel_mat4) <- i
    dfs <- rbind(dfs, sel_mat4)
    
    sel_mat3_af <- sel_combo[,2]
    sel_mat4_af <- data.frame(t(sel_mat3_af), check.names = FALSE)
    rownames(sel_mat4_af) <- i
    dfs_af <- rbind(dfs_af, sel_mat4_af)
    
    sel_mat3_dc <- sel_combo[,1]
    sel_mat4_dc <- data.frame(t(sel_mat3_dc), check.names = FALSE)
    rownames(sel_mat4_dc) <- i
    dfs_dc <- rbind(dfs_dc, sel_mat4_dc)
    
    
}

head(dfs)


# Order rows and make 1st row the colnames

colnames(dfs) <- dfs[1,]
dfs1 <- dfs[-1,]
dfs1$row <- as.integer(rownames(dfs1))
dfs1_ord <- dfs1[order(dfs1$row),]
cor <- dfs1_ord[,-c(dim(dfs1_ord)[2])]

colnames(dfs_af) <- dfs_af[1,]
dfsaf1 <- dfs_af[-1,]
dfsaf1$row <- as.integer(rownames(dfsaf1))
dfsaf1_ord <- dfsaf1[order(dfsaf1$row),]
afr <- dfsaf1_ord[,-c(dim(dfsaf1_ord)[2])]

colnames(dfs_dc) <- dfs_dc[1,]
dfsdc1 <- dfs_dc[-1,]
dfsdc1$row <- as.integer(rownames(dfsdc1))
dfsdc1_ord <- dfs1[order(dfsdc1$row),]
dcr <- dfsdc1_ord[,-c(dim(dfsdc1_ord)[2])]


head(cor)
head(afr)
head(dcr)

# Cell2mutation files inputted to DrugCell


write.csv(dfs, file = 'COMBO_vec.txt', quote = FALSE, row.names = FALSE)
write.csv(dfs_af, file = 'DEG_vec.txt', quote = FALSE, row.names = FALSE)
write.csv(dfs, file = 'MUT_vec.txt', quote = FALSE, row.names = FALSE)
# Used sed to remove header from csv files in Linux

# sed '1d' DEG_vec.txt > DEG_vec2.txt
# sed '1d' MUT_vec.txt > MUT_vec2.txt
# sed '1d' COMBO_vec.txt > COMBO_vec2.txt

#Vectors done

# filter out cell2ind.txt file (i.e some cell_lines in DrugCell and not within Affy Expression data)
dim(cor)
head(cell2ind)
colnames(cell2ind) <- c('dind', 'cell_line')
# create df that contains the cell_line indexs for the created vectors
ind <- data.frame(as.integer(rownames(dcr)))
colnames(ind) <- 'ind_vec'
dim(ind)

cell2ind_fil <- sqldf("SELECT * FROM ind INNER JOIN cell2ind ON ind.ind_vec = cell2ind.dind")

cell2ind_fil$ind <- c(0:914)

cell2ind_fil2 <- cell2ind_fil[,c(4,3)]

head(cell2ind_fil2)
dim(cell2ind_fil2)

library(readr)
# filtered cell2ind file that will be inputted to DrugCell
write_delim(cell2ind_fil2, file = 'cell2ind_v5.txt', delim = "\t", append = FALSE, col_names = FALSE, quote_escape = FALSE)



head(cell2ind_fil2)
dim(cell2ind_fil2)

# Load drugcell_all
drugcell_all <- read.delim(file ='drugcell_all.txt', check.names = FALSE, header = FALSE)
# Filter drugcell_all 
drugcell_all_fil <- sqldf("SELECT drugcell_all.V1, drugcell_all.V2, drugcell_all.V3 FROM cell2ind_fil2 INNER JOIN drugcell_all ON cell2ind_fil2.cell_line = drugcell_all.V1")

head(drugcell_all_fil)
dim(drugcell_all_fil)

# Load drugcell_test
drugcell_test <- read.delim(file ='drugcell_test.txt', check.names = FALSE, header = FALSE)
# Filter drugcell_test 
drugcell_test_fil <- sqldf("SELECT drugcell_test.V1, drugcell_test.V2, drugcell_test.V3 FROM cell2ind_fil2 INNER JOIN drugcell_test ON cell2ind_fil2.cell_line = drugcell_test.V1")

head(drugcell_test_fil)
dim(drugcell_test_fil)

# Load drugcell_train
drugcell_train <- read.delim(file ='drugcell_train.txt', check.names = FALSE, header = FALSE)
# Filter drugcell_train 
drugcell_train_fil <- sqldf("SELECT drugcell_train.V1, drugcell_train.V2, drugcell_train.V3 FROM cell2ind_fil2 INNER JOIN drugcell_train ON cell2ind_fil2.cell_line = drugcell_train.V1")

head(drugcell_train_fil)
dim(drugcell_train_fil)

# Load drugcell_train
drugcell_val <- read.delim(file ='drugcell_val.txt', check.names = FALSE, header = FALSE)
# Filter drugcell_train 
drugcell_val_fil <- sqldf("SELECT drugcell_val.V1, drugcell_val.V2, drugcell_val.V3 FROM cell2ind_fil2 INNER JOIN drugcell_val ON cell2ind_fil2.cell_line = drugcell_val.V1")

head(drugcell_val_fil)
dim(drugcell_val_fil)

write_delim(drugcell_val_fil, file = 'drugcell_val_fil_v5.txt', delim = "\t", append = FALSE, col_names = FALSE, quote_escape = FALSE)

dim(drugcell_val)

#write files out
write_delim(drugcell_test_fil, file = 'drugcell_test_fil_v5.txt', delim = "\t", append = FALSE, col_names = FALSE, quote_escape = FALSE)
write_delim(drugcell_train_fil, file = 'drugcell_train_fil_v5.txt', delim = "\t", append = FALSE, col_names = FALSE, quote_escape = FALSE)
write_delim(drugcell_all_fil, file = 'drugcell_all_fil_v5.txt', delim = "\t", append = FALSE, col_names = FALSE, quote_escape = FALSE)


# Double checking dim of testing files
library(readr)
#write_delim(cell2ind_fil2, file = 'cell2ind_v5.txt', delim = "\t", append = FALSE, col_names = FALSE, quote_escape = FALSE)

cell2ind_fil2 <- read.delim(file ='cell2ind_v5.txt', check.names = FALSE, header = FALSE)
head(cell2ind_fil2)
glimpse(cell2ind_fil2)
# Load drugcell_test
drugcell_test_new <- read.delim(file ='drugcell_test_new.txt', check.names = FALSE, header = FALSE)
dim(drugcell_test_new)
glimpse(drugcell_test_new)
# Filter drugcell_test 
drugcell_test_new_fil <- sqldf("SELECT drugcell_test_new.V1, drugcell_test_new.V2, drugcell_test_new.V3 FROM cell2ind_fil2 INNER JOIN drugcell_test_new ON cell2ind_fil2.V2 = drugcell_test_new.V1")
glimpse(drugcell_test_new_fil)
#head(drugcell_test_fil)
#dim(drugcell_test_fil)




