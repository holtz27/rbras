
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/rbras/main/aplica%C3%A7%C3%B5es/%5EBVSP.csv')
#ibovespa = read.csv('~/rstan/Aplicação/^BVSP.csv')
ibovespa = ibovespa[, c('Date', 'Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
View(ibovespa)
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
                                        real<lower=0.02,upper=0.5> lambda;
                                        real<lower=2> v; 
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
                                        // priori 1 
                                        //v ~ gamma(12, 0.8);
                                        // priori 2
                                        //target += 0.5 * log( v /( v+3 ) ) + 0.5 * log( trigamma( 0.5 * v ) - trigamma( 0.5 * (v+1) ) - 2 * ( v + 3 ) / ( v * square(v+1) ) );
                                        //priori 3
                                        lambda ~ uniform(0.02, 0.5);
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
                                          real<lower=0.02,upper=0.5> lambda;
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
                                          //priori 3
                                          lambda ~ uniform(0.02, 0.5);
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

N = 4e3
warmup = 1e4

fit0 = sampling(normal_fit, 
                list(T = length( log.ret ), 
                     y = log.ret,
                     y0 = 0),
                iter = warmup + N,
                warmup = warmup,
                chains = 4)
save(fit0, file = '~/rstan/Aplicação/fit0.RData')
load('~/rstan/Aplicação/fit0.RData')

parameter = extract( fit0 )































