
#######################################################################################################################


#Load packages

library("tidyverse")
library("dplyr")
library(data.table)
library(ggplot2)


##Data Preparation
#Read the data files (fang.mf and snp.mf) we are going to use.Because a header is present we will set the “header” argument to TRUE so that R knows the first line contains a header.

fang.mf <- read.delim("C:/Users/mfakude/OneDrive - Iowa State University/Documents/EEOB546_R_lesson/fang_et_al_genotypes.txt", header = TRUE)
snp.mf <- read.delim ("C:/Users/mfakude/OneDrive - Iowa State University/Documents/EEOB546_R_lesson/snp_position.txt", header = TRUE)

##Data inspection of fang_et_al_genotypes.txt (Now named fang.mf)

head(fang.mf) #shows first 6 rows
tail(fang.mf) #shows last 6 rows
dim(fang.mf) #returns the dimensions of data frame (number of rows and number of columns)
nrow(fang.mf) #number of rows
ncol(fang.mf) #number of columns
str(fang.mf) #structure of data frame - name, type and preview of data in each column
names(fang.mf) #shows the names attribute for a data frame, which gives the column names.
sapply(fang.mf, class) #shows the class of each column in the data frame
ls(fang.mf)#shows objects in the objects in the frame
View(fang.mf)#view the file in a spreadsheet-style

##Data inspection of snp_position.txt (Now named snp.mf)

head(snp.mf) #shows first 6 rows
tail(snp.mf) #shows last 6 rows
dim(snp.mf) #returns the dimensions of data frame (number of rows and number of columns)
nrow(snp.mf) #number of rows
ncol(snp.mf) #number of columns
str(snp.mf) #structure of data frame - name, type and preview of data in each column
names(snp.mf) #shows the names attribute for a data frame, which gives the column names
sapply(snp.mf, class) #shows the class of each column in the data frame
ls(fang.mf) #shows objects in the objects in the frame
View(snp.mf)
##Data Processing 
#Extract 3 columns (SNP_ID, Chromosome, Position) from snp.mf file and create a new file named snp_extracted

snp_extracted <- snp.mf[,c("SNP_ID", "Chromosome", "Position")]

###Save the newly created file containg the extracted columns and view it in a spreadsheet-style 

write.table(snp_extracted,"./Data_processing/snp_extracted.txt")
View(snp_extracted)

###Extract maize genotypes (ZMMIL, ZMMLR, ZMMMR) on 'Group' column from fang.mf file and create a new file named 'maize'.

maize <- fang.mf %>% filter(Group == "ZMMIL" | Group == "ZMMLR" | Group == "ZMMMR")

#Save the newly created file containg the extracted maize and view it

write.table(maize,"./Data_processing/maize.txt")
View(maize)
#Extract teosinte genotypes (ZMPBA, ZMPIL, ZMPJA) on 'Group' column from fang.mf file and create a new file named 'teosinte'.

teosinte <- fang.mf %>% filter(Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")

#Save the newly created file containg the extracted teosinte and view it

write.table(teosinte,"./Data_processing/teosinte.txt")
View(teosinte)

##Join object maize file and with snp_extracted file (For Maize) 
maize_transposed <- t(maize[,-c(1:3)]) %>% as.data.frame()
#Transpose maize file 
write.table(maize_transposed,"./Data_processing/maize_transposed.txt")
#Save the transposed maize file
maize_transposed$SNP_ID <- rownames(maize_transposed)
#
maize_final <- merge(snp_extracted,maize_transposed,by = 'SNP_ID')
#Join object maize file and with snp_extracted file
colnames(maize_final)[4:ncol(maize_final)] <- as.character(maize$Sample_ID)
#
maize_final[,-c(1:3)] <- lapply(maize_final[,-c(1:3)],as.character)
#
unique(maize_final$Chromosome)
#
write.table(maize_final,"./Data_processing/maize_final.txt")
#Save the joined maize file
View(maize_final)
#View the joined maize file

###Join object teosinte file and with snp_extracted file
#For teosinte
teosinte_transposed <- t(teosinte[,-c(1:3)]) %>% as.data.frame()
#Transpose maize file 
write.table(teosinte,"./Data_processing/teosinte_transposed.txt")
#Save the transposed maize file
teosinte_transposed$SNP_ID <- rownames(teosinte_transposed)
#
teosinte_final <- merge(snp_extracted,teosinte_transposed,by = 'SNP_ID')
#Join object teosinte file and with snp_extracted file  
colnames(teosinte_final)[4:ncol(teosinte_final)] <- as.character(teosinte$Sample_ID)
#
write.table(teosinte_final,"./Data_processing/teosinte_final.txt")
##Save the joined teosinte file

##Subseting maize chromosomes with SNPs ordered based on increasing position values and with missing data encoded by this symbol: ?
##files with chromosomes with SNPs arranged on increasing position values and with missing data encoded by this symbol '?' are name with an Q.
##files with chromosomes with SNPs arranged on decreasing position values and with missing data encoded by this symbol '-' are name with an H.
##For example, maize_chr1_Q or maize_chr1_H.

for (chr in 1:10) {
  maize_chr <- subset(maize_final,Chromosome == chr) %>% arrange(Position) 
  #subset by chromosome and arrange in ascending order
  maize_chr[maize_chr == '?/?'] <- '?' 
  #replace symbol ('?/?') with symbol ('?')
  write.table(maize_chr,file = paste("./Maize_chromosome/maize_chr",chr,"_Q.txt",sep = "")) 
  #save files and name each file by chromosome 
  maize_chr_d <- subset(maize_final,Chromosome == chr) %>% arrange(desc(Position))
  #subset by chromosome and arrange in descending order
  maize_chr_d[maize_chr_d == '?/?'] <- '-' 
  #replace symbol ('?/?') with symbol('?')
  write.table(maize_chr_d,file = paste("./Maize_chromosome/maize_chr",chr,"_H.txt",sep = ""))
  #save files and name each file by chromosome number and symbol  
  
  ##Subseting chromosomes for teosinte and arrange in ascending and decreasing order.
  teosinte_chr <- subset(teosinte_final,Chromosome == chr) %>% arrange(Position)
  #subset by chromosome and arrange in ascending order
  teosinte_chr[teosinte_chr == '?/?'] <- '?'
  #replace symbol ('?/?') with symbol ('?')
  write.table(teosinte_chr,file = paste("./Teosinte_chromosome/teosinte_chr",chr,"_Q.txt",sep = ""))
  #save files and name each file by chromosome 
  teosinte_chr_d <- subset(teosinte_final,Chromosome == chr) %>% arrange(desc(Position))
  #subset by chromosome and arrange in descending order
  teosinte_chr_d[teosinte_chr_d == '?/?'] <- '-'
  #replace symbol ('?/?') with symbol('?')
  write.table(teosinte_chr_d,file = paste("./Teosinte_chromosome/teosinte_chr",chr,"_H.txt",sep = ""))
  #save files and name each file by chromosome number and symbol
}
#1a.Plotting SNP counts per chromososome#
library()
snp.mf %>% 
  filter(Position < 10e1000) %>% 
  ggplot(aes(as.double(Chromosome))) +
  geom_bar(fill = 'purple', color = 'yellow') + 
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1) +
  scale_x_continuous(breaks = 1:10) +
  theme_replace() +
  ggtitle("snp count vs Chromosome") +
  ylab('Number of snp') +
  xlab('Chromosome') 
ggsave("./Data_visualization/Number of snp_chromosome.png")
#Filter the original SNP file by position and specify less than inf
#Pipe the filtering to create a bar graph of total counts of SNPs per chromosome
#Give the X-axis and the y-axis a label

#1b. Plotting SNP distribution across the chromosome#
snp.mf %>% filter(Position < 10e1000) %>% 
  ggplot(aes(as.double(Position)/1000000)) +
  geom_histogram(aes(y = ..density..), color = 'purple', fill = "purple", alpha = 0.4, bins = 20) + 
  geom_density(aes(as.double(Position)/1000000), color = "yellow") + 
  facet_wrap(~ as.double(Chromosome), scales = "free") +
  theme_replace() +
  ggtitle("SNP distribution across chromosomes") +
  xlab('Position (Mb)') +
  ylab('Density of SNPs')
ggsave("./Data_visualization/SNP_distribution_across_chromosomes.png")

#filter the position for anything less than infinity and map them to position/1000000, which is in double class
#Pipe this to create histogram SNP distribution in a chromosome and fill the histogram with green color
#Use facet Wrap to plot each chromosome
#Give the X-axis and the y-axis a label


##2.Missing data and amount of heterozygosity##

#2a. Proportion of homozygous and Hets by sample #
#Remove the 2nd tha the 3rd column in fang.mf and pipe it pivot_long by Sample_ID, to create a tibble of each SNP and its genotypic value 
#Mutate the file to add a colum with Homozygous and hetrozygous SNPs defined
mutate_Sampleid <- 
  fang.mf %>% select(-JG_OTU, -Group) %>%   
  pivot_longer(!Sample_ID) %>% 
  mutate(Allele = ifelse(value %in% c('C/C', 'G/G', 'A/A', 'T/T'), 'Homozygote', ifelse(value == '?/?', 'MD', 'Heterozygote')))  
View(mutate_Sampleid)

#Plotting by Sample#
library(RColorBrewer)
color_plots <- brewer.pal(3, "Set3")
mutate_genes %>% group_by(Sample_ID) %>%  count(Allele) %>% 
  ggplot(aes(fill = Allele, y = n, x = Sample_ID)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_plots) +
  ggtitle("Ratio of homozygotes, heterozygotes and missing data by genotype ") +
  ylab('Ratio') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("./Data_visualization/Gene_ratio by sample.png")

#2b. Proportion of homozygous and Hets by Group #
#Remove the 1st and the 3rd column in fang_et_al and pipe it to pivot_long by Group, to create a tibble of each group and its genotypic value 
#Mutate the file to add a column with Homozygous and Hetrozygous SNPs defined
mutate_Group <- 
  fang.mf %>% select(-JG_OTU, -Sample_ID) %>%   
  pivot_longer(!Group) %>% 
  mutate(Allele = ifelse(value %in% c('C/C', 'G/G', 'A/A', 'T/T'), 'Homozygote', ifelse(value == '?/?', 'MD', 'Heterozygote')))  
View(mutate_Group)

#### Plot by group ###
mutate_Group %>% group_by(Group) %>%  count(Allele) %>% 
  ggplot(aes(fill = Allele, y = n, x = Group)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color_plots) +
  ggtitle("Ratio of homozygotes, heterozygotes and missing data by Group ") +
  ylab('Ratio') 
ggsave("./Data_visualization/Ratio of zygosity and missing data.png")

##My_plot: Plotting by Sample ###

color2 <- brewer.pal(4, "Spectral")
mutate_groups %>% filter(Allele == "Homozygote") %>% group_by(Group) %>%  count(value) %>% 
  ggplot(aes(fill = value, y = n, x = Group)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = color2) +
  ggtitle("Ratio of homozygotes in each group") +
  ylab('Ratio') +
  theme_bw()
ggsave("./Data_visualization/my_plot.png")
