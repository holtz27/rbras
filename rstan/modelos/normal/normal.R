library(rstan)

options( mc.cores = 3 )
path = 'C:/Users/8936381/Documents/rstan/normal/stan_normal_model.stan'
model = stan_model( path )
################################################################################
# Running stan code
T = length( log.ret )
M  = 2e4
fit = sampling(model, 
               list(T = T, 
                    y = log.ret,
                    y0 = 0),
               iter = M,
               chains = 3)
#------------------ Save outputs
#saveRDS(parameters, file = 'C:/Users/8936381/Documents/rstan/normal/ibovespa.RData')
#parameters = readRDS('output_model_V')
#------------------ Plots
parameters = extract(fit)
x = summary(fit)
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2')
cols = c('mean', 'sd', 'X2.5.', 'X97.5.', 'n_eff','Rhat')
summary = Y[rows, cols]
summary$n_eff = 3 * M / (2 * summary$n_eff)
colnames(summary) = c('mean', 'sd', 'X2.5.', 'X97.5.', 'IF','Rhat')
summary
# h
rhat_h = Y[stringr::str_detect(row.names(Y), pattern = '^h'),][, 'Rhat']
plot(rhat_h, main = 'h')
abline(h = 0.95)
abline(h = 1)
abline(h = 1.05)
# Theta draws
draws = matrix(c( parameters$mu,
                  parameters$phi,
                  parameters$sigma,
                  parameters$b0, 
                  parameters$b1, 
                  parameters$b2 ), nrow = 6, byrow = TRUE)
draws = rbind( draws, t( as.matrix( parameters$h ) ) )
draws = rbind( draws, matrix(1, nrow = T, ncol = length(parameters$mu)) )
# plots
source('source/figures.R')
trace_plots(draws[1:6, ], names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2') )
abs_plots(draws[7:(T+6), ], draws[(T+7):nrow(draws), ], log.ret)
###############################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
###############################################################################
############### DIC
source('source/model_selection.R')
normal_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                        param_draws = draws[4:(nrow(draws)), ])
############### waic
normal_waic = svmsmn_waic(data = log.ret, y0 = 0,
                          draws = draws[4:(nrow(draws)), ])
############### loo
normal_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                        draws = draws[4:(nrow(draws)), ])
