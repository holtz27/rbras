# Reading data
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/rbras/main/aplica%C3%A7%C3%B5es/%5EBVSP.csv')
ibovespa = ibovespa[, c('Date', 'Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
#View(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) ) 
dates = as.Date( ibovespa[, 1] )

# Plots
library(ggplot2)
df = data.frame( Retorno = log.ret, Tempo = dates[-1] )

g = ggplot(df) + geom_line(aes(x = Tempo, y = Retorno))
g = g + scale_x_date(date_breaks = "58 month", date_labels = "%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size = 18),
                             axis.text.x = element_text(size = 16),
                             axis.text.y = element_text(size = 18))
g = g + xlab('')
g

h = ggplot( df, aes(Retorno) )
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x = element_text(size = 18),
                             axis.text.x = element_text(size = 18),
                             axis.text.y = element_text(size = 18))
h

gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2) 


data_summary = matrix(c( mean( log.ret ),
                         sd( log.ret ),
                         min( log.ret ),
                         max( log.ret ),
                         moments::skewness( log.ret ),
                         moments::kurtosis( log.ret ) ), nrow = 1)
colnames( data_summary ) = c( 'mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 2 )

Box.test(log.ret, lag = 12, type = 'Ljung-Box')
Box.test(log.ret**2, lag = 12, type = 'Ljung-Box')
acf(log.ret, lag.max = 5, plot = FALSE)


################################################################################
library(rstan)
options( mc.cores = 1 )
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
                                        //real<lower=0.025,upper=0.5> lambda;
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
                                        //lambda ~ uniform(0.025, 0.5);
                                        //v ~ exponential( lambda );
                                        v ~ gamma(2, 0.10);

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
                                          //real<lower=0.025,upper=2> lambda;
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
                                          //lambda ~ uniform(0.025, 2);
                                          //v ~ exponential( lambda );
                                          v ~ gamma(0.08, 0.04);
                                          
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
####################################
####################################
####################### normal
N = 4e3
warmup = 2e4
fit0 = sampling(normal_fit, 
                list(T = length( log.ret ), 
                     y = log.ret,
                     y0 = 0),
                iter = warmup + N,
                warmup = warmup,
                chains = 1)
#save(fit0, file = '~/rstan/Aplicação/fit0.RData')
load('~/rstan/Aplicação/fit0.RData')
parameters = extract( fit0 )
x = summary( fit0 )
Y = x$summary
#------------------ Tratando Y
#Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2')
cols = c('mean', 'sd', '2.5%', '97.5%')
summary = Y[rows, cols]
round( summary, digits = 3 )
# Theta draws
draws = matrix(c( parameters$mu,
                  parameters$phi,
                  parameters$sigma,
                  parameters$b0, 
                  parameters$b1, 
                  parameters$b2 ), nrow = 6, byrow = TRUE)
draws = rbind( draws, t( as.matrix( parameters$h ) ) )
draws = rbind( draws, matrix(1, nrow = T, ncol = length(parameters$mu)) )
############### Numeric Analysis
############################### theta
mcmc = coda::as.mcmc( t( draws[1:7, ] ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmc )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmc )
M / N_eff_theta
# plots
trace_plots(draws[1:6, ], names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2') )
abs_plots(draws[8:(T+6), ], log.ret)
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### dic
normal_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                        param_draws = draws[4:(nrow(draws)), ])
############### waic
normal_waic = svmsmn_waic(data = log.ret, y0 = 0,
                          draws = draws[4:(nrow(draws)), ])
############### loo
normal_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                        draws = draws[4:(nrow(draws)), ],
                        cores = 4)
####################################
####################################
####################### ts
N = 4e3
warmup = 2e4
fit1 = sampling(ts_fit, 
                list(T = length( log.ret ), 
                     y = log.ret,
                     y0 = 0),
                iter = warmup + N,
                warmup = warmup,
                chains = 1)
#save(fit0, file = '~/rstan/Aplicação/fit0.RData')
load('~/rstan/Aplicação/fit1.RData')
parameters = extract( fit1 )
x = summary(fit1)
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
cols = c('mean', 'sd', 'X2.5.', 'X50.','X97.5.')
summary = Y[rows, cols]
round( summary, digits = 3 )
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
############### Numeric Analysis
############################### theta
mcmc = coda::as.mcmc( t( draws[1:7, ] ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmc )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmc )
M / N_eff_theta
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
abs_plots(draws[9:(T+7), ],  log.ret)
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### dic
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
####################################
####################################
####################### slash
N = 4e3
warmup = 2e4
fit2 = sampling(slash_fit, 
                list(T = length( log.ret ), 
                     y = log.ret,
                     y0 = 0),
                iter = warmup + N,
                warmup = warmup,
                chains = 1)
#save(fit2, file = 'fit2.RData')
load('fit2.RData')
parameters = extract( fit2 )
x = summary( fit2 )
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
cols = c('mean', 'X2.5.', 'X97.5.')
summary = Y[rows, cols]
round( summary, digits = 3 )
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
#draws = draws[, seq(1, 4e4, 10)]
############### Numeric Analysis
############################### theta
mcmc = coda::as.mcmc( t( draws[1:7, ] ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmc )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmc )
N / N_eff_theta
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
abs_plots(draws[9:(T+7), ],  log.ret)
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### dic
slash_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                       param_draws = draws[ c( 4:6, 
                                               8:(T+7),
                                               (T+8):(nrow(draws)) ), ]
)
############### waic
slash_waic = svmsmn_waic(data = log.ret, y0 = 0,
                         draws = draws[ c( 4:6, 
                                           8:(T+7),
                                           (T+8):(nrow(draws)) ), ]
)
############### loo
slash_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                       draws = draws[ c( 4:6, 
                                         8:(T+7),
                                         (T+8):(nrow(draws)) ), ]
)
####################################
####################################
####################### vgamma
N = 4e3 
warmup = 2e4
fit3 = sampling(vgamma_fit, 
                list(T = length( log.ret ), 
                     y = log.ret,
                     y0 = 0),
                iter = warmup + N,
                warmup = warmup,
                chains = 1)
#save(fit3, file = '~/rstan/Aplicação/fit3.RData')
load('~/rstan/Aplicação/fit3.RData')
parameters = extract( fit3 )
x = summary( fit3 )
Y = x$summary
#------------------ Tratando Y
Y = data.frame(Y, row.names = row.names(Y))
#------------------ Avaliando convergência
rows = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
cols = c('mean', 'X2.5.', 'X97.5.')
summary = Y[rows, cols]
round( summary, digits = 3 )
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
############### Numeric Analysis
############################### theta
mcmc = coda::as.mcmc( t( draws[1:7, ] ) )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmc )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmc )
N_eff_theta
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
g = abs_plots(draws[8:(T+6), ], dates[-1], log.ret)
h = tail_plot(draws[(T+7):nrow(draws), ], dates[-1],'VG')
gridExtra::grid.arrange(h, g, nrow = 1, ncol = 2) 
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### dic
vgamma_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                        param_draws = draws[ c( 4:6, 
                                                8:(T+7),
                                                (T+8):(nrow(draws)) ), ]
)
############### waic
vgamma_waic = svmsmn_waic(data = log.ret, y0 = 0,
                          draws = draws[ c( 4:6, 
                                            8:(T+7),
                                            (T+8):(nrow(draws)) ), ]
)
############### loo
vgamma_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                        draws = draws[ c( 4:6, 
                                          8:(T+7),
                                          (T+8):(nrow(draws)) ), ]
)
