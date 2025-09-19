source(file = "./functions.R")

BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)

# get p
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

BRD$Number.of.Event.in.arm.2[BRD$Study.number == 2] <- BRD$Number.of.Event.in.arm.2[BRD$Study.number == 2] - 0.5

set.seed(20201221, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)

nrep <- 50000
n=c(50,100,150,200)
dat_para <- as.data.frame(n)
############## without previous information###########
for (c in 1:length(n)) {
  n_c <- n[c]
  
  pig_alloc_c <- rep(n_c/2,2)
  
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
    BRD_new_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_new_long, allstudies = T, sm = "OR")
    nma_sim_prev <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_pair,sm="OR",comb.fixed = T,comb.random = F)
    
    # lor_hat <- nma_sim_prev$TE.fixed[10,2]
    # se <- nma_sim_prev$seTE.fixed[10,2]
    # critical_vale <- abs(lor_hat)/se
    # pvalue <- 2*pnorm(-critical_vale)
    # if the simulate previous network shows the p-value is less than 0.1, then we conduct the new study, else not;
    if(nma_sim_prev$pval.fixed[10,2]<0.1){
      bioeq_single_s(p,p, pig_alloc_c)
    }
    
  }
  
  BE <- table_with_prev
  
  dat_para[c,2:3] <- pig_alloc_c
  dat_para[c,4:6] <- apply(BE, 2, sum)/nrep
  dat_para[c,7] <- nrow(BE)

}

write.csv(dat_para, file = "./res_error_single_dep.csv", row.names = F)
