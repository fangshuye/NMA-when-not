V3: modify the simulation process (also simulate the newwork)
    Only work with the effect size = 0

V4: simulate the whole network. Result of single study is copied from V3


V12: 
  - conduct the new study anyway: Res copy from V4
  - only conduct the new study when the previous network show that there are significant difference between the two treatment

V13:
  - conduct the new study anyway: Res copy from V4
  - only conduct the new study when the previous network show that there are significant difference between the two treatment: change a little code from V12 since the result of V12 is a little weird.


V15:
Recaculate the table of type I error, correct the calculation of SE
After checking, the SE is still big,
Test one of them, output all things.
res: dep

V16: 
check the SE again
Similar with V15, but got the result of table2: the new study will be conducted whatever the outcome of the existing network is
Res: indep

V17: 
change the comparison to TRIM and CEFTH.
P-value: 0.07731384 between two trts

V18: 
single (dep and indep)
Since the result in V17 is a little weried. Or we say, unstable. 10000 -> 50000

V19: 
same as V18. 50,000 -> 100,000

V20:
Add the result of 'conducting the new trial when the p-value of the previous results between 0.05 and 0.1'

