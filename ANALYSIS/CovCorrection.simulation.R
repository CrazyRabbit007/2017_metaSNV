n = 1000
res = 0 
div = runif(n, 0, 1)
mean.cov = 5 # Threshold ... 5 ?
sd.cov = 2 # Assumption : sd.cov = 2

for (i in 1:n){
  cov = rnorm(1, mean=mean.cov, sd=sd.cov)
  cov = max(cov,4) # SNVs not looked for if cov < 4
  res = res + div[i]*(cov/(cov-1))
}
res; print('Expected')
sum(div); print('Naive')
sum(div)*(mean.cov/(mean.cov-1)); print('Correction')
