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
                                        target += 0.5 * log( v /( v+3 ) ) + 0.5 * log( trigamma( 0.5 * v ) - trigamma( 0.5 * (v+1) ) - 2 * ( v + 3 ) / ( v * square(v+1) ) );
                                        // priori 3
                                        //lambda ~ uniform(0.02, 0.5);
                                        //v ~ exponential( lambda );
                                        
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
                                          // --- prioris
                                          target += 0.5 * log( v /( v+3 ) ) + 0.5 * log( trigamma( 0.5 * v ) - trigamma( 0.5 * (v+1) ) - 2 * ( v + 3 ) / ( v * square(v+1) ) );
                                          
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
                                          // --- prioris
                                          target += 0.5 * log( v /( v+3 ) ) + 0.5 * log( trigamma( 0.5 * v ) - trigamma( 0.5 * (v+1) ) - 2 * ( v + 3 ) / ( v * square(v+1) ) );
                                          
                                          //--- Sampling volatilitys:
                                          h_std ~ std_normal();
                                          for(t in 1:T) l[t] ~ inv_gamma( 0.5 * v, 0.5 * v );
                                          
                                          //--- Sampling observations:
                                          //y ~ normal( mu_t, exp(h/2) * l );
                                          for(t in 1:T) y[t] ~ normal( mu_t[t], exp( 0.5 * h[t] ) / sqrt( l[t] ) );
                                        }
                                        ' )
################################################################################
# Running stan code
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/svmsmn_stan_fitII.R' )
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/data/stable_data.R' )
source( 'https://raw.githubusercontent.com/holtz27/rbras/main/source/freq.models.R' )

df = list()
r = 5
set.seed( 753951 )
seeds = sample(1:1e6, r)
statistics = matrix(nrow = r, ncol = 6)
M = 7.5e3
alpha = 1.85
for(i in 1:r){
  if( i == 1 ) time.init = Sys.time()
  
  #data
  repeat{
    y = stable_data(mu = -1, phi = 0.985, sigma = 0.15,
                    b0 = 0.01, b1 = 0.1, b2 = -0.02,
                    y0 = 0,
                    a = alpha, # alpha e ( 0, 2 ] 
                    T = 2e3,
                    seed = seeds[ i ])
    if( moments::kurtosis(y) < 30 ) break
    else seeds[ i ] = sample( 1:1e6, 1 ) 
  }
  statistics[i, ] = matrix(c(mean(y),                   #Média
                             sd(y),                     #Desvio padrão
                             moments::skewness(y),      #Assimetria
                             moments::kurtosis(y),      #Curtose
                             min(y),                    #Mínimo
                             max(y) ),                   #Máximo
                           nrow = 1 )
  model_selection = data.frame()
  
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
  if( i == r ){
    time.final = Sys.time()
    colnames(statistics) = c('média', 'sd', 'skewnwss', 'kurtosis', 'min', 'max')
    time = time.final - time.init
    save(df, time, statistics, file = 'estudo2.RData')
  } 
  cat( '\r' )
}

#load('estudo2.RData')
statistics
freq.models( df )
