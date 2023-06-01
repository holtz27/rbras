library(rstan)

options( mc.cores = 4 )
path = 'C:/Users/8936381/Documents/rstan/Estudo_prioris/ts/ts_priori1.stan'
model = stan_model( path )
################################################################################
# Running stan code
source('C:/Users/8936381/Documents/source/ts_data.R')
T = 2e3
y = ts_data(mu = -1.0, phi = 0.985, sigma = 0.13,
            b0 = 0.01, b1 = 0.1, b2 = -0.23,
            y0 = 0,
            v = 10, 
            T = T,
            seed = 42
)
M  = 7500
fit = sampling(model, 
               list(T = T, 
                    y = y,
                    y0 = 0),
               iter = M,
               chains = 4 )
#------------------ Plots
parameters = extract(fit)
x = summary(fit)
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
cols = c('mean', 'sd', 'X2.5.', 'X97.5.', 'n_eff','Rhat')
summary = round( Y[rows, cols], digits = 3)
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
                  parameters$v ), nrow = 7, byrow = TRUE)
draws = rbind( draws, t( as.matrix( parameters$h ) ) )
draws = rbind( draws, t( as.matrix( parameters$l ) ) )
# jumps
lags = 5
jumps = seq(1, 0.5 * M * 4, by = lags)
draws = draws[, jumps]
################################################################################
############################## plots
source('C:/Users/8936381/Documents/source/figures.R')
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
abs_plots(draws[8:(T+7), ], draws[(T+8):nrow(draws), ],  y)
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### DIC
source('source/model_selection.R')
ts_dic = svmsmn_dic(data = y, y0 = 0, 
                    param_draws = draws[ c( 4:6, 
                                            8:(T+7),
                                            (T+8):(nrow(draws)) ), ]
                    )
############### waic
ts_waic = svmsmn_waic(data = y, y0 = 0,
                      draws = draws[ c( 4:6, 
                                        8:(T+7),
                                        (T+8):(nrow(draws)) ), ]
                      )
############### loo
ts_loo = svmsmn_loo(data = y, y0 = 0, 
                    draws = draws[ c( 4:6, 
                                      8:(T+7),
                                      (T+8):(nrow(draws)) ), ],
                    cores = 4
                    )






