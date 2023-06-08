library(rstan)
options( mc.cores = 4 )
################################################################################
# stan models
path0 = 'C:/Users/8936381/Documents/rstan/normal/stan_normal_model.stan'
normal_fit = stan_model( path0 )
# ts
path1 = 'C:/Users/8936381/Documents/rstan/ts/stan_ts_model.stan'
ts_fit = stan_model( path1 )
# slash
path2 = 'C:/Users/8936381/Documents/rstan/slash/stan_slash_model.stan'
slash_fit = stan_model( path2 )
# vgamma
path3 = 'C:/Users/8936381/Documents/rstan/vgama/stan_vgama_model.stan'
vgama_fit = stan_model( path3 )
################################################################################
# Running stan code
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/svmsmn_stan_fitII.R' )
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/stable_data.R' )

df = list()
r = 2
seeds = sample(1:1e6, r)
M = 7.5e3

for(i in 1:r){
  if( i == 1 ) time.init = Sys.time()
  #data
  y = stable_data(mu = -1, phi = 0.985, sigma = 0.13,
                  b0 = 0.001, b1 = 0.11, b2 = -0.21,
                  y0 = 0,
                  a = 1.85, # a e ( 0, 2 ] 
                  T = 2e3,
                  seed = seeds[ i ])
  
  model_selection = data.frame()
  #cat( '\r' )
  cat( paste0('réplica ', i ) )
  
  model_selection = svmsmn_stan_fitII(data = y,
                                      y0 = 0,
                                      model = normal_fit,
                                      M,
                                      nchains = 4, 
                                      lags = 5,
                                      model_name = 'normal',
                                      normal = TRUE,
                                      cores =  4)
  
  model_selection = rbind( model_selection, 
                           svmsmn_stan_fitII(data = y,
                                             y0 = 0,
                                             model = ts_fit,
                                             M,
                                             nchains = 4, 
                                             lags = 5,
                                             model_name = 'ts',
                                             normal = FALSE,
                                             cores =  4
                           )
  )
  
  model_selection = rbind( model_selection, 
                           svmsmn_stan_fitII(data = y,
                                             y0 = 0,
                                             model = slash_fit,
                                             M,
                                             nchains = 4, 
                                             lags = 5,
                                             model_name = 'slash',
                                             normal = FALSE,
                                             cores =  4)
  )
  
  model_selection = rbind( model_selection, 
                           svmsmn_stan_fitII(data = y,
                                             y0 = 0,
                                             model = vgama_fit,
                                             M,
                                             nchains = 4, 
                                             lags = 5,
                                             model_name = 'vgama',
                                             normal = FALSE,
                                             cores =  4)
  )
  
  df[[ i ]] = model_selection
  if( i == r ) time.final = Sys.time()
  cat( '\r' )
}

time.final - time.init
df
