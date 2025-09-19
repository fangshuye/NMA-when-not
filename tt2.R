information_size <- function(pc,pe){
  delta <- pc-pe
  p <- (pc+pe)/2
  sigma_2 <- p*(1-p)
  size <- 4*(qnorm(0.975)+qnorm(0.8))^2*sigma_2/(delta^2)
  round(size, 0)
}

alpha_function <- function(IF, alpha = 0.05){
  quantile <- qnorm(1-alpha/2)/sqrt(IF)
  2-2*pnorm(quantile)
}

information_size(0.3905,0.508)

# checked 
information_size(0.29,0.258) # the output is exactly the same as the TSA output
information_size(0.39,0.29)

information_size(0.3905,0.503)

# power.prop.test(power=0.8,p1=0.29,p2=0.258) # also simliar to this one # 6097 and 6099

# get alpha1
# param:
sample_size_exist <- 305+297
sample_size_new <- 100
p_control <- 0.30 # Gamithromycin
p_trt <- 0.25    # Florfenicol
IS <- information_size(0.29,0.258)
IF_1 <- sample_size_exist/IS
alpha_function(IF_1)
IF_2 <- (sample_size_exist+sample_size_new)/IS
alpha_function(IF_2) - alpha_function(IF_1)
obf.bd <- bounds(c(IF_1,IF_2),iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)

sample_size_exist <- 602
sample_size_new <- 103
IS <- information_size(0.39,0.29)
IF_1 <- sample_size_exist/IS
alpha_function(IF_1)
IF_2 <- (sample_size_exist+sample_size_new)/IS
alpha_function(IF_2) - alpha_function(IF_1)
obf.bd <- bounds(c(IF_1,IF_2),iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)

IS <- 1000
IF_1 <- sample_size_exist/IS
alpha_function(IF_1)
IF_2 <- (sample_size_exist+sample_size_new)/IS
alpha_function(IF_2) - alpha_function(IF_1)
obf.bd <- bounds(c(IF_1,IF_2),iuse=c(1,1),alpha=c(0.025,0.025))
summary(obf.bd)




sample_size_exist <- 305+297+100+100
sample_size_exist/IS
# confirm with the paper (discrete sequential boundaries for clinical trials)
alpha1 <- alpha_function(1/5)
c1 <- round(abs(qnorm(alpha1)),2) 
c1 # 4.23, same

alpha2 <- alpha_function(2/5) - alpha_function(1/5)
alpha2_adjust <- alpha2/(1-alpha1)
c2 <- round(abs(qnorm(alpha2_adjust)),2)
c2 # 2.89, same

alpha3 <- alpha_function(3/5) - alpha_function(2/5)  
alpha3_adjust <- alpha3/( pnorm(c1)*pnorm(c2) )
c3 <-  round(abs(qnorm(alpha3_adjust)),2)
c3
  
pnorm(2.41)*(1-pnorm(2.41))


