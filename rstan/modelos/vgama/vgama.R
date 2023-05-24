library(rstan)

options( mc.cores = 3 )
path = 'C:/Users/8936381/Documents/rstan/vgama/stan_vgama_model.stan'
model = stan_model( path )
################################################################################
# Running stan code
T = length( log.ret )
M  = 1e3
fit = sampling(model, 
               list(T = T, 
                    y = log.ret,
                    y0 = 0),
               iter = M,
               chains = 1 )
#------------------ Save outputs
#output = list( fit = fit, M = M )
#save(output, file = 'C:/Users/8936381/Documents/rstan/vgama/ibovespa_vgama.RData')
#load('C:/Users/8936381/Documents/rstan/ts/ibovespa_ts.RData')
#fit = output$fit
#M = output$M
#------------------ Plots
parameters = extract(fit)
x = summary(fit)
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
cols = c('mean', 'sd', 'X2.5.', 'X97.5.', 'n_eff','Rhat')
summary = Y[rows, cols]
summary
# h
rhat_h = Y[stringr::str_detect(row.names(Y), pattern = '^h') & 
           stringr::str_detect(row.names(Y), pattern = '_', negate = TRUE),][, 'Rhat']
plot(rhat_h, main = 'h')
abline(h = 0.95)
abline(h = 1)
abline(h = 1.05)
# l
rhat_l = Y[stringr::str_detect(row.names(Y), pattern = '^l'),][, 'Rhat']
plot(rhat_l, main = 'l')
abline(h = 0.95)
abline(h = 1)
abline(h = 1.05)
# Theta draws
draws = matrix(c( parameters$mu,
                  parameters$phi,
                  parameters$sigma,
                  parameters$b0, 
                  parameters$b1, 
                  parameters$b2,
                  parameters$v), nrow = 7, byrow = TRUE)
draws = rbind( draws, t( as.matrix( parameters$h ) ) )
draws = rbind( draws, t( as.matrix( parameters$l ) ) )
source('source/figures.R')
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
abs_plots(draws[8:(T+7), ], draws[(T+8):nrow(draws), ],  log.ret)
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### DIC
source('source/model_selection.R')
vgama_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                       param_draws = draws[ c( 4:6, 
                                            8:(T+7),
                                            (T+8):(nrow(draws)) ), ]
)
############### waic
vgama_waic = svmsmn_waic(data = log.ret, y0 = 0,
                         draws = draws[ c( 4:6, 
                                        8:(T+7),
                                        (T+8):(nrow(draws)) ), ]
)
############### loo
vgama_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                       draws = draws[ c( 4:6, 
                                      8:(T+7),
                                      (T+8):(nrow(draws)) ), ]
)



