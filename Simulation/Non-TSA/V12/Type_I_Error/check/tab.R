dat <- read.csv("check4.csv",header = T)

p_smller_0.1 <- dat[dat$V1<0.1,]
p_bigger_0.1 <- dat[dat$V1>=0.1,]

table(p_smller_0.1$V2)
table(p_bigger_0.1$V2)

apply(dat, 2, sum)/nrow(dat)
apply(dat, 2, mean)


apply(p_smller_0.1, 2, mean)

apply(p_bigger_0.1, 2, mean)
