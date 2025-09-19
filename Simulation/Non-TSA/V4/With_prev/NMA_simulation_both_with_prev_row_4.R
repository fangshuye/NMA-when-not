source(file = "./functions.R")

BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)



p_nac <- BRD_long %>% filter(t == "No active control") %>% summarise(p = sum(r)/sum(n)) %>% pull
lor_nac_2_TILD <- nma_old$TE.fixed[8, 10]
p <- lor2prob(p_nac,lor_nac_2_TILD)

# get the risk table
lor_nac_2_all <- nma_old$TE.fixed[8, ]
risk_all <- lor2prob(p_nac,lor_nac_2_all)
risk_all <- as.data.frame(risk_all)
risk_all$trt <- rownames(risk_all)
rownames(risk_all) <- rep(1:nrow(risk_all))
colnames(risk_all)[1] <- "p"

risk_all[risk_all$trt=="Ceftiofur pin","p"] <- risk_all[risk_all$trt=="Tildipirosin","p"]

nrep <- 10000

# need to change every time; from 1 to 4
r=4
n <- c(50,100,150,200)
n_r <- n[r]
dat_para <- as.data.frame(n_r)
############ with prev info##########################
pig_alloc_c <- rep(n_r/2,2)

BRD$Number.of.Event.in.arm.2[BRD$Study.number == 2] <- BRD$Number.of.Event.in.arm.2[BRD$Study.number == 2] - 0.5

set.seed(20201221, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)

table_with_prev <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  ### simulation the whole network###
  for (eachrow in 1:nrow(BRD)) {
    trt1 <- BRD$Arm.1[eachrow]
    trt2 <- BRD$Arm.2[eachrow]
    p1 <- risk_all[risk_all$trt==trt1,"p"]
    p2 <- risk_all[risk_all$trt==trt2,"p"]
    BRD[eachrow,4] <- rbinom(1, size =BRD[eachrow,7] , prob = p1)
    BRD[eachrow,5] <- rbinom(1, size =BRD[eachrow,8] , prob = p2)
    if(BRD$Number.of.arms[eachrow]==3){
      trt3 <- BRD$Arm.3[eachrow]
      p3 <- risk_all[risk_all$trt==trt3,"p"]
      BRD[eachrow,6] <- rbinom(1, size =BRD[eachrow,9] , prob = p3)
    }
  }
  BRD_new_long <- wide2long(BRD)
  
  bio_equal(p, p, sigma = 0, re_sigma = 0, pig_alloc_c, data_prev = BRD_new_long)
}

BE <- table_with_prev


dat_para[1,2:3] <- pig_alloc_c
dat_para[1,4:6] <- apply(BE, 2, sum)/nrep

write.csv(dat_para, file = paste0("./res_both_with_prev_row_",r,".csv"), row.names = F)


