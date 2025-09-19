library(network)
library(sna)
library(ggnetwork)
library(GGally)
library(patchwork)

generateMTMNetWork <- function(MTCdata, dataType = "Arm"){
  if( !(dataType %in% c("Arm", "Contrast") ) ){
    return(cat("Please input the correct data type, available options are \"Arm\" and \"Contrast\""))
  }
  maxArm <- max(MTCdata$Number.of.arms)
  if(dataType == "Arm"){
    ans <- data.frame(matrix(ncol = 4, nrow = 0))
    names(ans) <- c("t1", "t2", "n1", "n2")
    for(i in 2:maxArm){
      tmp_all <- MTCdata[MTCdata$Number.of.arms == i, ]
      for(j in 1:(i-1)){
        for(k in (j+1):i){
          tmp <- tmp_all[, c(paste0("Arm", j), paste0("Arm", k), paste0("Total.number.in.arm.", j), 
                             paste0("Total.number.in.arm.", k))]
          names(tmp) <- c("t1", "t2", "n1", "n2")
          ans <- rbind(ans, tmp)
        }
      }
    }
    ## order the data 
    ans <- ans[order(ans$n1), ]
  }
  if(dataType == "Contrast"){
    ans <- data.frame(matrix(ncol = 2, nrow = 0))
    names(ans) <- c("t1", "t2")
    for(i in 2:maxArm){
      tmp_all <- MTCdata[MTCdata$Number.of.arms == i, ]
      for(j in 1:(i-1)){
        for(k in (j+1):i){
          tmp <- tmp_all[, c(paste0("Arm", j), paste0("Arm", k))]
          names(tmp) <- c("t1", "t2")
          ans <- rbind(ans, tmp)
        }
      }
    }
  }
  rownames(ans) <- NULL
  return(ans)
}

mtmNetWorkNodeName <- function(MTCdata, map_txt, arm_name = "abbr", save_res = F){
  if(!(arm_name %in% c("abbr", "full"))){return(cat("Please input correct arm name type, available options are \"abbr\" and \"full\""))}
  if(arm_name == "abbr"){
    treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(abbr) %>% dplyr::pull()
  }
  if(arm_name == "full"){
    treatNms <- map_txt %>% dplyr::arrange(number) %>% dplyr::select(name) %>% dplyr::pull()
  }
  ## no. of treatment 
  NT <- length(treatNms)
  ## treatment frequency - the col numbers are the cols with the numerical txt indictors 
  maxArm <- max(MTCdata$Number.of.arms) ## get the maximum number of arms in a study
  ntFreq <- NULL
  for(i in 1:NT){
    ntFreq[i] <- (MTCdata[,c(paste0("Arm", 1:maxArm))] == i) %>% sum(na.rm = T)
  }
  ## combine treatment name with frequency 
  treat_new <- ""
  for(i in 1:NT){
    treat_new[i] <- paste(treatNms[i], "(", ntFreq[i], ")", sep = "")
  }
  
  ## write treatment frequency out 
  num_trea <- cbind(treatNms,ntFreq)
  if(save_res == T){
    write.csv(num_trea, file = "./data/treatmentFrequencySummary.csv", row.names = FALSE)
  }
  numTrials <- data.frame(matrix(nrow = 1, ncol = maxArm-1))
  for(i in 2:maxArm){
    numTrials[1,(i-1)] <- sum(MTCdata$Number.of.arms == i)
  }
  colnames(numTrials) <- paste0("num", 2:maxArm, "ArmsTrials")
  
  return(list(treat_new = treat_new, 
              treat_freq = ntFreq,
              totalTrialArms  = sum(ntFreq),
              numTrials = numTrials))
}

mat_treat.fun <- function(t1,t2,percomparison,nameoftreatments,VAR1,VAR2,graphtitle,thickness,nodetextsize,nodesize,
                          vertex.col = "red",...) {
  if (percomparison==T){#if the first two arguments are the treatments and the 
    #second two arguments are the corresponding randomized sample sizes
    t1<-t1
    t2<-t2
  }
  #in this part we transform the two arguments: ID and Treatment
  #in the same as previous form
  else{
    ID<-t1
    Treat<-t2    
    treat <- as.double(as.factor(Treat))
    id <- ID
    nt <- length(unique(treat))
    ns <- length(unique(ID))
    na <- table(match(id, unique(id)))
    #_____Check whether the treatments are 1 to nt
    if(max(sort(unique(treat)) - c(1:nt)) > 0) {
      treat <- as.factor(treat)
      levels(treat) <- c(1:nt)
      treat <- as.double(treat)
    }
    TT <- matrix(NA, nrow = ns, ncol = max(na))
    
    for(i in 1:ns) {
      TT[i, 1:na[i]] <- treat[id == unique(id)[i]]
    }
    u <- TT
    new.id <- 1:length(table(id))  
    TT <- (apply(TT, 1, sort))
    na <- na
    u <- c()
    torepeat <- na[na > 2]
    # if there are more than 3 arms
    if(sum(na > 2))
    {
      for(i in 1:(sum(na > 2))) 
      {
        u <- rbind(u, t(combn(unlist(TT[na > 2][i]),2)))
      }
    }
    
    TT<- rbind(matrix(unlist(TT[na == 2]), ncol = 2, byrow = T), u)
    t1<-TT[,1]
    t2<-TT[,2]    
  }# end of the else that transforms (ID,Treatment) to (t1,t2)
  #_____________________________________________________________________
  
  numoftreatments=max(cbind(t1,t2));
  
  # Default adjacency-style matrix initialization
  mat_treat<-matrix(0,numoftreatments,numoftreatments);
  
  if (missing(nameoftreatments))
  {
    nam="";
    for(i in 1:numoftreatments) 
    { 
      nam[i] <- paste("treat",i,sep=".")
      assign(nam[i],1:i)
    }
    colnames(mat_treat)<-nam;
    rownames(mat_treat)<-nam;
  }
  else
  {
    colnames(mat_treat)<-nameoftreatments;
    rownames(mat_treat)<-nameoftreatments;
  }
  
  # Adjusted adjacency-style matrix initialization (for thickness of links)
  thickness_new<-matrix(0,numoftreatments,numoftreatments);
  
  if (missing(nameoftreatments))
  {
    nam="";
    for(i in 1:numoftreatments) 
    { 
      nam[i] <- paste("treat",i,sep=".")
      assign(nam[i],1:i)
    }
    colnames(thickness_new)<-nam;
    rownames(thickness_new)<-nam;
  }
  else
  {
    colnames(thickness_new)<-nameoftreatments;
    rownames(thickness_new)<-nameoftreatments;
  }
  #_____________________________________________________________________
  # based on the frequency of treatments we construct the node's thickness
  # divisor has fixed value equal to 5 for aesthetical reasons
  if (missing(nodesize)){
    divisor=5} else {
      divisor=5/nodesize
    }
  
  if (missing(VAR1)){
    nodethickness<-table(c(t1,t2))/divisor }else{
      nodethickness<-VAR1
    }
  
  if (missing(thickness)){
    thickness=10
  }
  
  #combined treatments (for default link thickness)
  tr<-cbind(t1,t2)
  
  for (i in 1:length(t1))
  {
    mat_treat[tr[i,1],tr[i,2]]<-mat_treat[tr[i,1],tr[i,2]]+1
  }
  
  mat_treatb<-t(mat_treat)
  mat_treat<-mat_treat+mat_treatb
  
  return(mat_treat)
}


setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/data')
MTCdata <- read.csv(file = "updated_dataset.csv", stringsAsFactors = F)

map_txt <- openxlsx::read.xlsx(xlsxFile = "./mapping.xlsx")

pairwise_comp <- generateMTMNetWork(MTCdata, dataType = "Arm")
networkNode <-  mtmNetWorkNodeName(MTCdata, map_txt = map_txt)
treat_new <- networkNode$treat_new



mat_treat <- mat_treat.fun(pairwise_comp$t1, pairwise_comp$t2, percomparison = T, 
                           nameoftreatments = treat_new)
mat_treat_tmp <- ifelse(mat_treat > 0 , 0.5, 0)
# TRIM and CEFTH
# need modification if the comparison changes.
mat_treat_tmp[12,2] <- 2
mat_treat_tmp[2,12] <- 2


### add network attributes
tmp <- network(mat_treat_tmp)
set.edge.value(tmp, "thick", mat_treat_tmp)
#ggnet2(tmp, mode = 'circle', label = T, edge.size = "thick")
# name the 2th and 12th trt
network.vertex.names(tmp) <- c(rep(NA, 1),
                               'CEFTH',
                               rep(NA,9),
                               "TRIM", 
                               rep(NA, 1))

tmp %v% "color" = c(rep('grey', 1),
                    'black',
                    rep('grey',9),
                    "black", 
                    rep('grey', 1))

### adjust nodes positions
x = gplot.layout.circle(tmp, NULL)


# switch TRMI and CEFTP (12 and 3)
# make TRMI to the right side
TRMI_p = x[12,]
x[12,] = x[3,]
x[3,] = TRMI_p


for (i in 1:13) {
  if(i!=12 & i!=2){
    x[i,] <- 0.75*x[i,]
  }
}

tmp %v% 'x' = x[,1]
tmp %v% 'y' = x[,2]


tmp %v% "t1" = c(rep(1,13),NA)

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/Network_plot/')
##### plot the network
ggnet2(tmp, mode = c("x", "y"), color = 'color', label = F, edge.color = c("color", "grey"), edge.size = "thick")  +
  guides(size = FALSE)  +
  annotate("text", x = c(0.846,1), y = c(1,0.848), label = c("CEFTH",'TRIM')) +
  # legend line: grey
  annotate("segment", x = 0, xend = 0.12, y = 1, yend = 1,  colour = "black", size = 2)+
  # legend line: black
  annotate("segment", x = 0, xend = 0.12, y = 0.95, yend = 0.95,  colour = "grey", size = 0.5)+
  # legend text
  annotate("text", x = c(0.21, 0.257), y = c(1, 0.95), label = c("The new trial", "The existing network")) 

# ggsave(paste0("BMC_figure_V2",".png"),
#        width = 6.5,  # 2.63 - 7.5
#        height = 6.5, # <8.75
#        dpi = 300,
#        units = "in")
ggsave(paste0("figure1",".jpeg"),
       width = 180,  # 2.63 - 7.5
       height = 180, # <8.75
       dpi = 300,
       units = "mm")
#############################################
################# no new trial ##############
#############################################

mat_treat <- mat_treat.fun(pairwise_comp$t1, pairwise_comp$t2, percomparison = T, 
                           nameoftreatments = treat_new)
mat_treat_tmp <- ifelse(mat_treat > 0 , 0.5, 0)
# TRIM and CEFTH
# need modification if the comparison changes.
mat_treat_tmp[12,2] <- 2
mat_treat_tmp[2,12] <- 2


### add network attributes
tmp <- network(mat_treat_tmp)
set.edge.value(tmp, "thick", mat_treat_tmp)
#ggnet2(tmp, mode = 'circle', label = T, edge.size = "thick")
# name the 2th and 12th trt
network.vertex.names(tmp) <- c(rep(NA, 1),
                               'CEFTH',
                               rep(NA,9),
                               "TRIM", 
                               rep(NA, 1))
tmp %v% "color" = rep('grey', 13)
# tmp %v% "color" = c(rep('grey', 1),
#                     'black',
#                     rep('grey',9),
#                     "black", 
#                     rep('grey', 1))

### adjust nodes positions
x = gplot.layout.circle(tmp, NULL)


# switch TRMI and CEFTP (12 and 3)
# make TRMI to the right side
TRMI_p = x[12,]
x[12,] = x[3,]
x[3,] = TRMI_p


for (i in 1:13) {
  if(i!=12 & i!=2){
    x[i,] <- 0.75*x[i,]
  }
}

tmp %v% 'x' = x[,1]
tmp %v% 'y' = x[,2]


tmp %v% "t1" = c(rep(1,13),NA)

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/Network_plot/')
##### plot the network
ggnet2(tmp, mode = c("x", "y"), color = 'color', label = F, edge.color = c("color", "grey"), edge.size = "thick")  +
  guides(size = FALSE) 
ggsave(paste0("BMC_figure_no_label",".png"),
       width = 6.5,  # 2.63 - 7.5
       height = 6.5, # <8.75
       dpi = 300,
       units = "in")


ggnet2(tmp, mode = c("x", "y"), color = 'color', label = F, edge.color = c("color", "grey"), edge.size = "thick")  +
  guides(size = FALSE)  +
  annotate("text", x = c(0.8,0.95), y = c(1,0.84), label = c("CEFTH",'TRIM'))

ggsave(paste0("BMC_figure_V2_no_new_sim1",".png"),
       width = 6.5,  # 2.63 - 7.5
       height = 6.5, # <8.75
       dpi = 300,
       units = "in")

#############################################
################# sim2 ##############
#############################################

mat_treat <- mat_treat.fun(pairwise_comp$t1, pairwise_comp$t2, percomparison = T, 
                           nameoftreatments = treat_new)
mat_treat_tmp <- ifelse(mat_treat > 0 , 0.5, 0)
# GAMI anf FLOR
# need modification if the comparison changes.
# mat_treat_tmp[7,8] <- 2
# mat_treat_tmp[8,7] <- 2


### add network attributes
tmp <- network(mat_treat_tmp)
set.edge.value(tmp, "thick", mat_treat_tmp)
#ggnet2(tmp, mode = 'circle', label = T, edge.size = "thick")
# name the 2th and 12th trt
network.vertex.names(tmp) <- c(rep(NA, 6),
                               'FLOR',
                               "GAMI", 
                               rep(NA, 5))
tmp %v% "color" = rep('grey', 13)
# tmp %v% "color" = c(rep('grey', 1),
#                     'black',
#                     rep('grey',9),
#                     "black", 
#                     rep('grey', 1))

### adjust nodes positions
x = gplot.layout.circle(tmp, NULL)

tmp %v% 'x' = x[,1]
tmp %v% 'y' = x[,2]


tmp %v% "t1" = c(rep(1,13),NA)

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/Network_plot/')
##### plot the network
ggnet2(tmp, mode = c("x", "y"), color = 'color', label = F, edge.color = c("color", "grey"), edge.size = "thick")  +
  guides(size = FALSE)  +
  annotate("text", x = c(0.3,0.7), y = c(0.007,0.007), label = c("GAMI",'FLOR'))

ggsave(paste0("BMC_figure_V2_no_new_sim2",".png"),
       width = 6.5,  # 2.63 - 7.5
       height = 6.5, # <8.75
       dpi = 300,
       units = "in")

#############################################
################# sim2 with new trial ##############
#############################################

mat_treat <- mat_treat.fun(pairwise_comp$t1, pairwise_comp$t2, percomparison = T, 
                           nameoftreatments = treat_new)
mat_treat_tmp <- ifelse(mat_treat > 0 , 0.5, 0)
# GAMI anf FLOR
# need modification if the comparison changes.
mat_treat_tmp[7,8] <- 2
mat_treat_tmp[8,7] <- 2


### add network attributes
tmp <- network(mat_treat_tmp)
set.edge.value(tmp, "thick", mat_treat_tmp)
#ggnet2(tmp, mode = 'circle', label = T, edge.size = "thick")
# name the 2th and 12th trt
network.vertex.names(tmp) <- c(rep(NA, 6),
                               'FLOR',
                               "GAMI", 
                               rep(NA, 5))
tmp %v% "color" = c(rep('grey', 6),
                    'black',"black",
                    rep('grey',5))

### adjust nodes positions
x = gplot.layout.circle(tmp, NULL)


tmp %v% 'x' = x[,1]
tmp %v% 'y' = x[,2]


tmp %v% "t1" = c(rep(1,13),NA)

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/Network_plot/')
##### plot the network
ggnet2(tmp, mode = c("x", "y"), color = 'color', label = F, edge.color = c("color", "grey"), edge.size = "thick")  +
  guides(size = FALSE)  +
  annotate("text", x = c(0.3,0.7), y = c(0.007,0.007), label = c("GAMI",'FLOR'))

ggsave(paste0("BMC_figure_V2_new_sim2",".png"),
       width = 6.5,  # 2.63 - 7.5
       height = 6.5, # <8.75
       dpi = 300,
       units = "in")

##### plot the network
x = gplot.layout.circle(tmp, NULL)
x[13,2] <- 0.755

tmp %v% 'x' = x[,1]
tmp %v% 'y' = x[,2]


tmp %v% "t1" = c(rep(1,13),NA)

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/Network_plot/')
##### plot the network
ggnet2(tmp, mode = c("x", "y"), color = 'color', label = F, edge.color = c("color", "grey"), edge.size = "thick")  +
  guides(size = FALSE)  +
  # add grey line between GAMI and FLOR
  annotate("segment", x = 0.403, xend = 0.6, y = 0.01, yend = 0.01,  colour = "grey", size = 0.5)+
  annotate("text", x = c(0.3,0.7), y = c(0.007,0.007), label = c("GAMI",'FLOR'))+
  # legend line: grey
  annotate("segment", x = 0, xend = 0.12, y = 1, yend = 1,  colour = "black", size = 2)+
  # legend line: black
  annotate("segment", x = 0, xend = 0.12, y = 0.95, yend = 0.95,  colour = "grey", size = 0.5)+
  # legend text
  annotate("text", x = c(0.21, 0.257), y = c(1, 0.95), label = c("The new trial", "The existing network")) 

# ggsave(paste0("BMC_figure_sim2",".png"),
#        width = 6.5,  # 2.63 - 7.5
#        height = 6.5, # <8.75
#        dpi = 300,
#        units = "in")

ggsave(paste0("figure2",".jpeg"),
       width = 180,  # 2.63 - 7.5
       height = 180, # <8.75
       dpi = 300,
       units = "mm")

