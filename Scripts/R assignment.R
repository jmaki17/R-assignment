##Start of R assignment

##I started out by creating a new project in R linking it to my GitHub in the shell
##I then moved the fang et al genotypes file and the snp_position files into the R project folder 
fang <- read.table("fang_et_al_genotypes.txt", fill=TRUE, header = TRUE, sep="\t", stringsAsFactors = FALSE)
snps <- read.table("snp_position.txt", fill=TRUE, header=TRUE, sep="\t", stringsAsFactors = FALSE)
##loading the fang and snp files into my global environment
library(ggplot2)
library(dplyr)
##loading up the packages needed for this assignment
str(fang)
str(snps)
##checking the structure of "fang" and "snps"
head(fang)
head(snps)
##checking the header of both fang and snps
typeof(fang)
typeof(snps)
##determining the type of both files ()
length(fang)
length(snps)
##determining the length of both files
class(fang)
class(snps)
##confirming that both files are being recognized as dataframes
dim(fang)
dim(snps)
##determining the number of columns and rows in each of the files
names(fang)
names(snps)
##determining the column names to 
maize<-filter(fang, fang$Group=="ZMMIL"| fang$Group=="ZMMLR" | fang$Group=="ZMMMR")
teosinte<-filter(fang, fang$Group=="ZMPBA"| fang$Group=="ZMPIL" | fang$Group=="ZMPJA")
##extracting the maize and teosinte genotypes from the fang dataset using filter
View(maize)
View(teosinte)
##view the results to inspect the output and insure it is right
trans_maize <- t(maize)
trans_teosinte <- t(teosinte)
##transposing the maize and teosinte, makes the dataframe into a matrix
df_trans_maize <- as.data.frame(trans_maize, stringsAsFactors = FALSE)
df_trans_teosinte <- as.data.frame(trans_teosinte, stringsAsFactors = FALSE)
##converting the resulting matrix into a dataframe
merged_maize <- merge(snps, df_trans_maize, by.x=1, by.y=0)
merged_teosinte <- merge(snps, df_trans_teosinte, by.x=1, by.y=0)
##merging the snps and transposed maize and teosinte data frames
View(merged_maize)
View(merged_teosinte)
##viewing my merged files to ensure they look correct
ordered_merged_maize <- merged_maize[,c(1,3,4,2,5:ncol(merged_maize))]
ordered_merged_teosinte <- merged_teosinte[,c(1,3,4,2,5:ncol(merged_teosinte))]
##ordering the columns so that the first column is SNP_ID, the second is Chromosome
##and the third column is Position and all the genotype data is behind that
##probably a bit fragile of a solution (need to know the number of variables in each data frame)
View(ordered_merged_maize)
View(ordered_merged_teosinte)
##viewing ordered dataframes to ensure they look correct
hyph_maize <-lapply(ordered_merged_maize, gsub, pattern = "?", replacement = "-", fixed = TRUE)
hyph_teosinte <-lapply(ordered_merged_teosinte, gsub, pattern = "?", replacement = "-", fixed = TRUE)
##beginning the substitution of "-" in the place of "?" in both data frames, but this causes
##the data frame to get converted to a value

df_hyph_maize <- data.frame(hyph_maize, stringsAsFactors = FALSE)
df_hyph_teosinte <- data.frame(hyph_teosinte, stringsAsFactors = FALSE)
##converting the maize and teosinte "-" values back into data frames
for (i in 1:10){
  ordered_merged_maize_sub=subset(ordered_merged_maize,ordered_merged_maize$Chromosome==i)
  maize_chromos<- ordered_merged_maize_sub[order(strtoi(ordered_merged_maize_sub$Position)),]
write.table(maize_chromos, file =sprintf("maize_chromo%i.txt",i), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)}
subset
##creating text files for the maize genotypes with ? for unknown values in ascending order
##This creates the files in the working directory for the project (can be seen in git tab to right)
for (i in 1:10){
hyph_maize_sub =subset(df_hyph_maize, df_hyph_maize$Chromosome==i)
  hyph_maize_chromos<- hyph_maize_sub[order(strtoi(hyph_maize_sub$Position), decreasing = TRUE),]
  write.table(hyph_maize_chromos, file =sprintf("hyph_maize_chromo%i.txt",i), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)}
##creating text files for the maize genotypes with "-" for unknown values in descending order
for (i in 1:10){
  ordered_merged_teosinte_sub=subset(ordered_merged_teosinte,ordered_merged_teosinte$Chromosome==i)
  teosinte_chromos<- ordered_merged_teosinte_sub[order(strtoi(ordered_merged_teosinte_sub$Position)),]
  write.table(teosinte_chromos, file =sprintf("teosinte_chromo%i.txt",i), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)}

##creating text files for the teosinte genotypes with ? for unknown values in ascending order
for (i in 1:10){
  hyph_teosinte_sub =subset(df_hyph_teosinte, df_hyph_teosinte$Chromosome==i)
  hyph_teosinte_chromos<- hyph_teosinte_sub[order(strtoi(hyph_teosinte_sub$Position), decreasing = TRUE),]
  write.table(hyph_teosinte_chromos, file =sprintf("hyph_teosinte_chromo%i.txt",i), sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)}
##creating text files for the teosinte genotypes with "-" for unknown values in descending order

##Assignment part 2 electric boogaloo

##graph #1
maize_and_teosinte_groups <-filter(fang, fang$Group=="ZMMIL"| fang$Group=="ZMMLR" | fang$Group=="ZMMMR" | fang$Group=="ZMPBA"| fang$Group=="ZMPIL" | fang$Group=="ZMPJA")
##pulling out all of the group data for both maize and teosinte groups
df_m_and_t_groups <- as.data.frame(t(maize_and_teosinte_groups), stringsAsFactors = FALSE)
##transposing the maize and teosinte groups and making them into a data frame
merged_m_and_t <- merge(snps, df_m_and_t_groups, by.x=1, by.y=0)
##merging the maize and teosinte groups with the snps dataset
install.packages("reshape2")
library(reshape2)
##installing and loading the "reshape2" package
##melting the merged and ordered maize and teosinte files to make life easier for ggplot
melt_merged_m_and_t <- melt(merged_m_and_t)

df_m_and_t <-as.data.frame(table(melt_merged_m_and_t$Chromosome))
##making the melted and merged file into a dataframe, generates a dataframe which gives
##the chromosome number as well as the frequency of snps on each chromosome
##this will give us the total number of SNPs in our dataset on each chromosome
ggplot(df_m_and_t)+geom_col(aes(x=df_m_and_t$Var1, y=df_m_and_t$Freq))+ xlab("Chromosome") +
  ylab("Number of SNPs")
##actually creating the plot of snps/chromosome
##it appears as though chromosome 1 accounts for the most SNPs within the dataset for maize and teosinte

##determining number of SNPs from each group
snps_per_group <-data.frame(table(maize_and_teosinte_groups$Group))
View(snps_per_group)
##to determine the number of snps that are contributed by specific groups, the combined
##maize and teosinte groups file was converted to a dataframe that gives the number 
##of SNPs contributed to the total by each group
##veiwing the group, we see that ZMMLR contributes 1256 of the SNPs, making it the largest contributor of SNPs to the dataset

##Plot of my choice
ggplot(snps_per_group) + geom_col(aes(x=snps_per_group$Var1, y=snps_per_group$Freq)) +
  xlab("Group") + ylab("Number of SNPs") 
##made a bar graph that shows the number of SNPs contributed by each group

##Plot #2
##I am reading this question to mean we're looking at the all the snips at 1 position
##of our choosing (like abph1.20 for instance)
plot2_m_and_t <- maize_and_teosinte_groups[,-2]
##getting rid of the second column from the maize and teosinte dataset because it's not needed
melt_plot2 <- melt(plot2_m_and_t,id=c("Sample_ID","Group"))
##melting the dataset
colnames(melt_plot2)[3:4] <- c("snp","snp_nucleotide") 
##changing the names of the second and third columns to remember more easily
my_fun <- function(x) {
  if ( x == "A/A" | x == "T/T" | x == "C/C" | x == "G/G") {
    return(1)
  }
  else if (x == "?/?") {
    return(0)
  }
  else {return(2)}
}
##creating a function that we can use in conjunction with lapply to sort our newly created
##"snp call id

melt_plot2$new <- lapply(melt_plot2$snp_nucleotide,my_fun)
##creating a "new" column which contains information on whether the sample is homozygous(1),
##heterozygous(2), or an unknown value (0)
View(melt_plot2)
##view the plot to make sure it worked
melt_plot2$new <- as.numeric(melt_plot2$new)
##coverting the data in the "new" column into numerics
arr_melt_plot2 <- melt_plot2 %>% arrange(Group, Sample_ID)
##arranging the dataframe based firstly on group, then by sample id
plot2_group_table <- arr_melt_plot2 %>% group_by(arr_melt_plot2$Group,arr_melt_plot2$Sample_ID,arr_melt_plot2$new) %>% summarize(n=n())
##calculating the number of snps within each group
plot2_total_table <- arr_melt_plot2 %>% group_by(arr_melt_plot2$Group,arr_melt_plot2$Sample_ID) %>% summarize(total = n())
##calculating the total number of SNPs in the dataset to allow us to determine percentages
merged_plot2 <- merge(plot2_group_table,plot2_total_table,by.x=2, by.y=2)
##merging the 2 dataframes using the sample id column
View(merged_plot2)
##checking the dataframe to make sure it worked
merged_plot2$Percent <- merged_plot2$n/merged_plot2$total
##calculated the percentage of each snp compared to the total numbers of snps and 
##added them to a new row in the table called "Percent"

ggplot(data = merged_plot2, aes(x = merged_plot2$`arr_melt_plot2$Sample_ID`, y = merged_plot2$Percent, colour = merged_plot2$`arr_melt_plot2$Group.x`, xlab("Sample ID", ylab("Percentage")))) +
  legend(title = "Group") + facet_grid(~ merged_plot2$`arr_melt_plot2$new`)
##generating the plot for the graph

ggplot(arr_melt_plot2, aes(x = arr_melt_plot2$new)) +geom_bar(aes(fill=arr_melt_plot2$Group)) + 
  xlab("Heterozygosity (0 = NA, 1 = Homozygous, 2 = Heterozygous)") + 
  ylab("SNPs")
