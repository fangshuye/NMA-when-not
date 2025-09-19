library(ldbounds)

t <- seq(0.2,1,length=5)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)
plot(obf.bd)


t <- c(0.0987,0.1151)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)


t <- c(0.2404153,0.3202875)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)


t <- c(0.2404153,0.2803514)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)

t <- c(602/702,1)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)


t <- c(602/1000,702/1000)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd) # 2.66， 2.49

t <- c(602/800,702/800)
obf.bd <- bounds(t,iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd) # 2.34， 2.20




