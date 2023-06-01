library(rstan)
options( mc.cores = 4 )
source('C:/Users/8936381/Documents/source/estudo_priori_stan.R')
# n° replicações
r = 50
seeds = sample(1:1e6, r)
# priori 1
path = 'C:/Users/8936381/Documents/rstan/Estudo_prioris/ts/ts_priori1.stan'
model1 = stan_model( path )
priori1 = estudo_priori_stan(y0 = 0,
                             model = model1,
                             M = 2.5e3,
                             nchains = 4,
                             lags = 1,
                             r = r,
                             seeds = seeds )
priori1$summary
priori1$time
