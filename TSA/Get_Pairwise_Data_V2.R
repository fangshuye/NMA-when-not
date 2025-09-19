setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/NMA/')
source(file = "./functions_when_not.R")

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/')
BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
#BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
#BRD_long <- wide2long(BRD_r)
#BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
#nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)

#p_value <- nma_old$pval.fixed
#p_value[p_value < 0.05] <- '1'
# write.csv(p_value,
#           'p_value_nma.csv',
#           row.names = T)
library(dplyr)
BRD_nac <- BRD %>% filter(Arm2 == 1)

BRD_sub <- BRD %>% 
  filter(Arm.1 == 'Gamithromycin') %>%
  filter(Arm.2 == 'Tilmicosin')

BRD_sub <- BRD %>% 
  filter(Arm.1 == 'Tilmicosin') %>%
  filter(Arm.2 == 'Gamithromycin')




pair_tab <- BRD %>%
  group_by(Arm.1,Arm.2) %>%
  summarize(Study_num = n())
table(BRD_nac$Arm.1)

pair_tab <- pair_tab %>%
  arrange(desc(Study_num))
trt1 <- pair_tab$Arm.1[1] # ENFO
trt2 <- pair_tab$Arm.2[2] # NAC
# they have 21 pairwise comparisons

BRD_pair <- BRD_r %>% 
  filter(Arm.1 == trt1) %>%
  filter(Arm.2 == trt2)

BRD_pair$Study.number <- rownames(BRD_pair)
wide_long_pair <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$Study.number,N)
  t <- NULL
  n <- NULL
  r <- NULL
  for(i in c(2,1)){
    r <- c(r, eval(parse(text = paste0("MTCdata$Number.of.Event.in.arm.",i, sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$Total.number.in.arm.",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
  }
  res <- data.frame(id = study, t = t, r = r, n = n)
  res <- res %>% dplyr::filter(!is.na(n)) %>% arrange(id)
  res
}

BRD_pair_long <- wide_long_pair(BRD_pair)

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/data')
write.csv(BRD_pair, file = paste0("./pairwise_dataset",".csv"), row.names = F)
write.csv(BRD_pair_long, file = paste0("./pairwise_long_dataset",".csv"), row.names = F)

# convert to TSA format
rm(list = ls())
library(dplyr)
library(tibble)
library(stringr)
setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/TSA')
samplecol_dat <- read.csv('RevMan_Output.csv',
                          header = T,
                          check.names = F,
                          na.strings=" ")

samplecol_dat$`Group Label 1`[2] <- 'Enrofloxacin'
samplecol_dat$`Group Label 2`[2] <- 'No active control'
samplecol_dat$Name[1:2] <- c('ENFO vs NAC','BRD')
samplecol_dat <- samplecol_dat[c(1,2),]
# write.csv(test, 
#           file = paste0("./RevMan_Output_copy",".csv"), 
#           na = '',
#           row.names = F,
#           quote = F)
# version 1 
setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/data')
dat <- read.csv('pairwise_dataset.csv',header = T)
colnames(dat)[2] <- 'Year of study'
colnames(dat)[3] <- 'Name'
colnames(dat)[4:5] <- c('Events 1','Events 2')
colnames(dat)[7:8] <- c('Total 1','Total 2')

dat <- dat[,c(2:5,7:8)]
dat[is.na(dat)] <- rep(1999,4)
dat$Name <- str_replace_all(dat$Name, ",", "-")
# add the new study data
num_row <- nrow(dat)
dat[num_row+1,] <- c(2021,NA,10,46,90,50)
dat$Name[num_row+1] <- 'New Study'

dat <- dat %>%
  add_column('SD 1' = 0,
             'Mean 1' = 0,
             'SD 2' = 0,
             'Mean 2' = 0,
             'Data Type' = ' ',
             'Group Label 1' = ' ',
             'Group Label 2' = ' ',
             'Comparison Number' = 1,
             'Outcome Number' = 1,
             'Subgroup Number' = 0)

colname_need <- colnames(samplecol_dat)
dat <- dat[,colname_need]

combine_dat <- rbind(samplecol_dat,dat)
combine_dat$`Events 1`[2] <- sum(dat$`Events 1`)
combine_dat$`Events 2`[2] <- sum(dat$`Events 2`)
combine_dat$`Total 1`[2] <- sum(dat$`Total 1`)
combine_dat$`Total 2`[2] <- sum(dat$`Total 2`)


write.csv(combine_dat, 
          file = paste0("./pairwise_dataset_TSA_format",".csv"), 
          na = '',
          row.names = F,
          quote = F)


# version 2 
setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/data')
dat <- read.csv('pairwise_dataset_simulated.csv',header = T)
colnames(dat)[2] <- 'Year of study'
colnames(dat)[3] <- 'Name'
colnames(dat)[4:5] <- c('Events 1','Events 2')
colnames(dat)[6:7] <- c('Total 1','Total 2')

dat[is.na(dat)] <- rep(1999,4)
dat$Name <- str_replace_all(dat$Name, ",", "-")

dat <- dat %>%
  add_column('SD 1' = 0,
             'Mean 1' = 0,
             'SD 2' = 0,
             'Mean 2' = 0,
             'Data Type' = ' ',
             'Group Label 1' = ' ',
             'Group Label 2' = ' ',
             'Comparison Number' = 1,
             'Outcome Number' = 1,
             'Subgroup Number' = 0)

colname_need <- colnames(samplecol_dat)
dat <- dat[,colname_need]

combine_dat <- rbind(samplecol_dat,dat)
combine_dat$`Events 1`[2] <- sum(dat$`Events 1`)
combine_dat$`Events 2`[2] <- sum(dat$`Events 2`)
combine_dat$`Total 1`[2] <- sum(dat$`Total 1`)
combine_dat$`Total 2`[2] <- sum(dat$`Total 2`)


write.csv(combine_dat, 
          file = paste0("./pairwise_dataset_simulated_TSA_format",".csv"), 
          na = '',
          row.names = F,
          quote = F)










