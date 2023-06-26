
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/rbras/main/aplica%C3%A7%C3%B5es/%5EBVSP.csv')
#ibovespa = read.csv('~/rstan/Aplicação/^BVSP.csv')
ibovespa = ibovespa[, c('Date', 'Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
ibovespa = tail(ibovespa, n = 2000)
#View(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) ) 
dates = as.Date( ibovespa[, 1] )

plot( log.ret ~ dates[-1], type = 'l' )

summary = matrix(c( mean( log.ret ),
                    sd( log.ret ),
                    min( log.ret ),
                    max( log.ret ),
                    moments::skewness( log.ret ),
                    moments::kurtosis( log.ret ) ), nrow = 1)
colnames( summary ) = c( 'mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
summary

Box.test(log.ret, lag = 12, type = 'Ljung-Box')
Box.test(log.ret**2, lag = 12, type = 'Ljung-Box')
acf(log.ret, lag.max = 5, plot = FALSE)
################################################################################
library(rstan)
options( mc.cores = 4 )
################################################################################
# stan models
normal_fit = stan_model( model_code = 'data{
                                              int<lower=0> T;
                                              real y0;
                                              vector[T] y;
                                            }
                                            parameters{
                                              real mu;                         // mean log volatility
                                              real<lower=-1,upper=1> phiT;     // persistence of volatility
                                              real<lower=0> s2;
                                              vector[T] h_std;                 // std log volatility time t
                                              real<lower=-1,upper=1> b1T;
                                              real b0;
                                              real b2;
                                            }
                                            transformed parameters{
                                              real<lower=-1,upper=1> phi;
                                              real<lower=-1,upper=1> b1;
                                              real<lower=0> sigma;
                                              vector[T] h;  // now h ~ normal(0, sigma)
                                              vector[T] mu_t;
                                              phi = (2*phiT - 1);
                                              b1 = (2*b1T - 1);
                                              sigma = sqrt(s2);
                                              
                                              //--- Volatilitys:
                                              h = h_std * sigma;
                                              h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
                                              h += mu;
                                              for (t in 2:T) {
                                                h[t] += phi * (h[t - 1] - mu);
                                              }
                                              
                                              //--- Means:
                                              mu_t[1] = b0 + y0*b1 + exp( h[1] ) * b2;
                                              for(t in 2:T){
                                                mu_t[t] = b0 + b1 * y[t-1]  + exp( h[t] ) * b2;
                                              }
                                            }
                                            model{
                                              // --- prioris
                                              mu ~ normal(0, 3.2);
                                              phiT ~ beta(20, 1.5);
                                              b0 ~ normal(0, 3.2);
                                              b1T ~ beta(5, 1.5);
                                              b2 ~ normal(0, 3.2);
                                              s2 ~ inv_gamma(2.5, 0.025);
                                              
                                              //--- Sampling volatilitys:
                                              h_std ~ std_normal();
                                              
                                              //--- Sampling observations:
                                              y ~ normal( mu_t, exp(h/2) );
                                            }
                                            ' )

# ts
ts_fit = stan_model( model_code = 'data{
                                        int<lower=0> T;
                                        real y0;
                                        vector[T] y;
                                      }
                                      parameters{
                                        real mu;                         // mean log volatility
                                        real<lower=-1,upper=1> phiT;     // persistence of volatility
                                        real<lower=0> s2;
                                        vector[T] h_std;  // std log volatility time t
                                        vector<lower=0>[T] l;
                                        real<lower=0.025,upper=0.5> lambda;
                                        real<lower=2, upper=40> v; 
                                        real<lower=-1,upper=1> b1T;
                                        real b0;
                                        real b2;
                                      }
                                      transformed parameters{
                                        real<lower=-1,upper=1> phi;
                                        real<lower=-1,upper=1> b1;
                                        real<lower=0> sigma;
                                        vector[T] h;  // now h ~ normal(0, sigma)
                                        vector[T] mu_t;
                                        phi = (2*phiT - 1);
                                        b1 = (2*b1T - 1);
                                        sigma = sqrt(s2);
                                        
                                        //--- Volatilitys:
                                        h = h_std * sigma;
                                        h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
                                        h += mu;
                                        for (t in 2:T) {
                                          h[t] += phi * (h[t - 1] - mu);
                                        }
                                        
                                        //--- Means:
                                        mu_t[1] = b0 + y0 * b1 + exp( h[1] ) * b2;
                                        for(t in 2:T){
                                          mu_t[t] = b0 + b1 * y[t-1]  + exp( h[t] ) * b2;
                                        }
                                      }
                                      model{
                                        // --- prioris
                                        mu ~ normal(0, 3.2);
                                        phiT ~ beta(20, 1.5);
                                        b0 ~ normal(0, 3.2);
                                        b1T ~ beta(5, 1.5);
                                        b2 ~ normal(0, 3.2);
                                        s2 ~ inv_gamma(2.5, 0.025);
                                        //priori 3
                                        lambda ~ uniform(0.025, 0.5);
                                        v ~ exponential( lambda );
                                        
                                        //--- Sampling volatilitys:
                                        h_std ~ std_normal();
                                        for(t in 1:T) l[t] ~ gamma(0.5 * v, 0.5 * v);
                                        
                                        //--- Sampling observations:
                                        //y ~ normal( mu_t, exp(h/2) * l );
                                        for(t in 1:T) y[t] ~ normal( mu_t[t], exp( 0.5 * h[t] ) / sqrt( l[t] ) );
                                      }
                                      ' )
# slash
slash_fit = stan_model( model_code = 'data{
                                          int<lower=0> T;
                                          real y0;
                                          vector[T] y;
                                        }
                                        parameters{
                                          real mu;                         // mean log volatility
                                          real<lower=-1,upper=1> phiT;     // persistence of volatility
                                          real<lower=0> s2;
                                          vector[T] h_std;  // std log volatility time t
                                          vector<lower=0, upper=1> [T] l;
                                          real<lower=0> v; 
                                          real<lower=-1,upper=1> b1T;
                                          real b0;
                                          real b2;
                                        }
                                        transformed parameters{
                                          real<lower=-1,upper=1> phi;
                                          real<lower=-1,upper=1> b1;
                                          real<lower=0> sigma;
                                          vector[T] h;  // now h ~ normal(0, sigma)
                                          vector[T] mu_t;
                                          phi = (2*phiT - 1);
                                          b1 = (2*b1T - 1);
                                          sigma = sqrt(s2);
                                          
                                          //--- Volatilitys:
                                          h = h_std * sigma;
                                          h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
                                          h += mu;
                                          for (t in 2:T) {
                                            h[t] += phi * (h[t - 1] - mu);
                                          }
                                          
                                          //--- Means:
                                          mu_t[1] = b0 + y0 * b1 + exp( h[1] ) * b2;
                                          for(t in 2:T){
                                            mu_t[t] = b0 + b1 * y[t-1]  + exp( h[t] ) * b2;
                                          }
                                        }
                                        model{
                                          // --- prioris
                                          mu ~ normal(0, 3.2);
                                          phiT ~ beta(20, 1.5);
                                          b0 ~ normal(0, 3.2);
                                          b1T ~ beta(5, 1.5);
                                          b2 ~ normal(0, 3.2);
                                          s2 ~ inv_gamma(2.5, 0.025);
                                          //priori 1 
                                          v ~ gamma(0.08, 0.04);
                                          
                                          //--- Sampling volatilitys:
                                          h_std ~ std_normal();
                                          for(t in 1:T) l[t] ~ beta( v, 1 );
                                          
                                          //--- Sampling observations:
                                          //y ~ normal( mu_t, exp(h/2) * l );
                                          for(t in 1:T) y[t] ~ normal( mu_t[t], exp( 0.5 * h[t] ) / sqrt( l[t] ) );
                                        }
                                        ' )
# vgamma
vgama_fit = stan_model( model_code = 'data{
                                          int<lower=0> T;
                                          real y0;
                                          vector[T] y;
                                        }
                                        parameters{
                                          real mu;                         // mean log volatility
                                          real<lower=-1,upper=1> phiT;     // persistence of volatility
                                          real<lower=0> s2;
                                          vector[T] h_std;  // std log volatility time t
                                          vector<lower=0> [T] l;
                                          real<lower=0.025,upper=2> lambda;
                                          real<lower=0, upper=40> v; 
                                          real<lower=-1,upper=1> b1T;
                                          real b0;
                                          real b2;
                                        }
                                        transformed parameters{
                                          real<lower=-1,upper=1> phi;
                                          real<lower=-1,upper=1> b1;
                                          real<lower=0> sigma;
                                          vector[T] h;  // now h ~ normal(0, sigma)
                                          vector[T] mu_t;
                                          phi = (2*phiT - 1);
                                          b1 = (2*b1T - 1);
                                          sigma = sqrt(s2);
                                          
                                          //--- Volatilitys:
                                          h = h_std * sigma;
                                          h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
                                          h += mu;
                                          for (t in 2:T) {
                                            h[t] += phi * (h[t - 1] - mu);
                                          }
                                          
                                          //--- Means:
                                          mu_t[1] = b0 + y0 * b1 + exp( h[1] ) * b2;
                                          for(t in 2:T){
                                            mu_t[t] = b0 + b1 * y[t-1]  + exp( h[t] ) * b2;
                                          }
                                        }
                                        model{
                                          // --- prioris
                                          mu ~ normal(0, 3.2);
                                          phiT ~ beta(20, 1.5);
                                          b0 ~ normal(0, 3.2);
                                          b1T ~ beta(5, 1.5);
                                          b2 ~ normal(0, 3.2);
                                          s2 ~ inv_gamma(2.5, 0.025);
                                          //priori 3
                                          lambda ~ uniform(0.025, 2);
                                          v ~ exponential( lambda );
                                          
                                          //--- Sampling volatilitys:
                                          h_std ~ std_normal();
                                          for(t in 1:T) l[t] ~ inv_gamma( 0.5 * v, 0.5 * v );
                                          
                                          //--- Sampling observations:
                                          //y ~ normal( mu_t, exp(h/2) * l );
                                          for(t in 1:T) y[t] ~ normal( mu_t[t], exp( 0.5 * h[t] ) / sqrt( l[t] ) );
                                        }
                                        ' )
################################################################################
source('https://raw.githubusercontent.com/holtz27/rbras/main/source/figures.R')
source('https://raw.githubusercontent.com/holtz27/rbras/main/source/model_selection.R')
N = 4e3 / 4
warmup = 2e4
#fit0 = sampling(normal_fit, 
#                list(T = length( log.ret ), 
#                     y = log.ret,
#                     y0 = 0),
#                iter = warmup + N,
#                warmup = warmup,
#                chains = 4)
#save(fit0, file = '~/rstan/Aplicação/fit0.RData')
#load('~/rstan/Aplicação/fit0.RData')
#parameters = extract( fit0 )
#x = summary(fit0)
#Y = x$summary
#------------------ Tratando Y
#Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
#rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2')
#cols = c('mean', 'sd', 'X2.5.', 'X97.5.', 'n_eff','Rhat')
#summary = Y[rows, cols]
#summary
# h
#rhat_h = Y[stringr::str_detect(row.names(Y), pattern = '^h') & 
#             stringr::str_detect(row.names(Y), pattern = '_', negate = TRUE),][, 'Rhat']
#plot(rhat_h, main = 'h')
#abline(h = 0.95)
#abline(h = 1)
#abline(h = 1.05)
# Theta draws
#draws = matrix(c( parameters$mu,
#                  parameters$phi,
#                  parameters$sigma,
#                  parameters$b0, 
#                  parameters$b1, 
#                  parameters$b2 ), nrow = 6, byrow = TRUE)
#draws = rbind( draws, t( as.matrix( parameters$h ) ) )
#draws = rbind( draws, matrix(1, nrow = T, ncol = length(parameters$mu)) )
# plots
#trace_plots(draws[1:6, ], names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2') )
#abs_plots(draws[7:(T+6), ], log.ret)
# model selection
#normal_dic = svmsmn_dic(data = log.ret, y0 = 0, 
#                        param_draws = draws[4:(nrow(draws)), ])
############### waic
#normal_waic = svmsmn_waic(data = log.ret, y0 = 0,
#                          draws = draws[4:(nrow(draws)), ])
############### loo
#normal_loo = svmsmn_loo(data = log.ret, y0 = 0, 
#                        draws = draws[4:(nrow(draws)), ],
#                        cores = 4)
####################################
####################################
####################### ts
fit1 = sampling(ts_fit, 
                list(T = length( log.ret ), 
                     y = log.ret,
                     y0 = 0),
                iter = warmup + N,
                warmup = warmup,
                chains = 4)
#save(fit0, file = '~/rstan/Aplicação/fit0.RData')
#load('~/rstan/Aplicação/fit0.RData')
parameters = extract( fit1 )
x = summary(fit1)
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
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
abs_plots(draws[8:(T+7), ],  log.ret)
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### DIC
ts_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                    param_draws = draws[ c( 4:6, 
                                            8:(T+7),
                                            (T+8):(nrow(draws)) ), ]
)
############### waic
ts_waic = svmsmn_waic(data = log.ret, y0 = 0,
                      draws = draws[ c( 4:6, 
                                        8:(T+7),
                                        (T+8):(nrow(draws)) ), ]
)
############### loo
ts_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                    draws = draws[ c( 4:6, 
                                      8:(T+7),
                                      (T+8):(nrow(draws)) ), ]
)





























