library(rstan)
options( mc.cores = 4 )
source('https://raw.githubusercontent.com/holtz27/rbras/main/source/estudo_priori_stan.R')
# priori 1
path = 'C:/Users/8936381/Documents/rstan/ts/priori2.stan'
model = stan_model( path )
# n° replicações
r = 1
seeds = sample(1:1e6, r)
priori2 = estudo_priori_stan(y0 = 0,
                             model = model,
                             M = 7.5e3,
                             nchains = 4,
                             lags = 5,
                             r = r,
                             seeds = seeds )
priori1$summary
priori1$time
