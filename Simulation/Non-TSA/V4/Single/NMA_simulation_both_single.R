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


cl <- makeCluster(16)
registerDoParallel(cl)
set.seed(20200927)

nrep <- 10000
n=c(50,100,150,200)
dat_para <- as.data.frame(n)
############## without previous information###########
for (c in 1:length(n)) {
  n_c <- n[c]
  
  pig_alloc_c <- rep(n_c/2,2)
  
  table_wo_prev <- foreach(i = 1:nrep) %dorng%{
    bioeq_single_s(p,p, pig_alloc_c)
  }
  
  BE <- do.call("rbind", table_wo_prev)
  
  dat_para[c,2:3] <- pig_alloc_c
  dat_para[c,4:6] <- apply(BE, 2, sum)/nrep

}
stopCluster(cl)

write.csv(dat_para, file = "./res_both_single.csv", row.names = F)
