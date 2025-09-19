Infor <- function(a,b,c,d){
  1/(1/a+1/b+1/c+1/d)
}
LIL <- function(z_unadjust, cumu_infor, lambda = 2){
  deno <- sqrt(lambda*log(log(cumu_infor)))
  no <- z_unadjust
  round(no/deno,2)
}

setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_II_When_not/data/V2/Simulated')
dat <- read.csv('pairwise_dataset_simulated_115.csv')
colnames(dat)[1] <- 'Study.number'
dat$'Arm.1' <- 'Florfenicol'
dat$'Arm.2' <- 'Gamithromycin'
setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/NMA/')
source(file = "./functions_when_not.R")

wide_long_pair <- function(MTCdata){
  N <- 2
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

# 1st and 2nd study
dat_long <- wide_long_pair(dat)
dat_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = dat_long, allstudies = T, sm = "OR")
nma_2 <- netmeta(TE,seTE,treat1,treat2,studlab,data=dat_pair,sm="OR",comb.fixed = T,comb.random = F)
lor_hat <- nma_2$TE.fixed[2,1]
se <- nma_2$seTE.fixed[2,1]
z <- lor_hat/se
round(z,2) # checked
# -0.91
LIL(z, 1/se^2, 2) # under checking: correct: -1.03

# 1st study
dat_long <- wide_long_pair(dat[1,])
dat_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = dat_long, allstudies = T, sm = "OR")
nma_2 <- netmeta(TE,seTE,treat1,treat2,studlab,data=dat_pair,sm="OR",comb.fixed = T,comb.random = F)
#nma_2 <- netmeta(TE,seTE,treat1,treat2,studlab,data=dat_pair,sm="RD",comb.fixed = T,comb.random = F)

lor_hat <- nma_2$TE.fixed[2,1]
se <- nma_2$seTE.fixed[2,1]
z <- lor_hat/se
round(z,2) # checked
info1 <- 1/se^2

#info <- (qnorm(0.975)+qnorm(0.8))^2/lor_hat^2
# -1.3
LIL(z, info1, 2) # under checking: correct: -1.44
LIL(z, info, 2) # under checking: correct: -1.44

# 2nd study
dat_long <- wide_long_pair(dat[2,])
dat_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = dat_long, allstudies = T, sm = "OR")
nma_2 <- netmeta(TE,seTE,treat1,treat2,studlab,data=dat_pair,sm="OR",comb.fixed = T,comb.random = F)
lor_hat <- nma_2$TE.fixed[2,1]
se <- nma_2$seTE.fixed[2,1]
z <- lor_hat/se
round(z,2) # checked
info2 <- 1/se^2


# checked, the same
#Infor(75,297-75,100,305-100)
#1/se^2


