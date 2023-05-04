library(ggplot2)
library(rstan)
library(stringr)

options( mc.cores = 3 )
#getwd()
path = 'C:/Users/8936381/Documents/rstan/normal/stan_normal_model.stan'
model = stan_model( path )
################################################################################
#### Lendo os dados...

y = log.ret[1:3000]
par( mfrow = c(1, 2) )
plot(y, type = 'l')
hist(y, breaks = 40)
par( mfrow = c(1, 1) )
#Descritiva dos log-retornos
mean(y)                  #Média
median(y)                #Mediana
sd(y)                    #Desvio padrão
moments::skewness(y)     #Assimetria
moments::kurtosis(y)     #Curtose
min(y)                   #Mínimo
max(y)                   #Máximo

# Running stan code
T = length(y)
M  = 5e3
fit = sampling(model, 
               list(T = T, 
                    y = y,
                    y0 = 0),
               iter = M,
               chains = 3)

#------------------ Save outputs
#saveRDS(parameters, file = 'output_model_V')
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
data = Y[rows, cols]
data



# h
rhat_h = Y[str_detect(row.names(Y), pattern = '^h'),][, 'Rhat']
plot(rhat_h, main = 'h')
abline(h = 0.95)
abline(h = 1)
abline(h = 1.05)

# mu
mu = parameters$mu
par(mfrow=c(1,3))
plot(mu, type='l', main='', xlab='')
hist(mu, main = round(mean(mu), 3), xlab='', col = 'white', breaks = 40)
abline(v = mean(mu), col = 'blue', lwd=2, lty=2)
abline(v = Y['mu', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['mu', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(mu, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1))  

# phi
phi = parameters$phi
par(mfrow=c(1,3))
plot(phi, type='l', main='', xlab='')
hist(phi, main = round(mean(phi), 3), xlab='', col = 'white', breaks = 40)
abline(v = mean(phi), col = 'blue', lwd=2, lty=2)
abline(v = Y['phi', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['phi', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(phi, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1))

# sigma
s = parameters$sigma
par(mfrow=c(1,3))
plot(s, type='l', main='', xlab='')
hist(s, main = round(mean(s), 3), xlab='', col = 'white', breaks = 40)
abline(v = mean(s), col = 'blue', lwd=2, lty=2)
abline(v = Y['sigma', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['sigma', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(s, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 

################################################################################
############### Análise numérica
theta = matrix( c( mu, phi, s ), ncol = 3 )
mcmcchain_theta = coda::as.mcmc( theta )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_theta = coda::geweke.diag( mcmcchain_theta )
CD_theta
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_theta = coda::effectiveSize( mcmcchain_theta )
IF_theta = nrow(theta) / N_eff_theta
IF_theta
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_theta = round( apply( theta, MARGIN = 2, FUN = sd) / sqrt( N_eff_theta ), 5 )
mc_error_theta
################################################################################

# b0
b0 = parameters$b0
par(mfrow=c(1,3))
plot(b0, type='l', main='', xlab='')
hist(b0, main=round(mean(b0), 8), xlab='', col = 'white', breaks = 40)
abline(v = mean(b0), col = 'blue', lwd=2, lty=2)
abline(v = Y['b0', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['b0', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(b0, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 

# b1
b1 = parameters$b1
par(mfrow=c(1,3))
plot(b1, type='l', main='', xlab='')
hist(b1, main=round(mean(b1), 8), xlab='', col = 'white', breaks = 40)
abline(v = mean(b1), col = 'blue', lwd=2, lty=2)
abline(v = Y['b1', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['b1', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(b1, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 

# b2
b2 = parameters$b2
par(mfrow=c(1,3))
plot(b2, type='l', main='', xlab='')
hist(b2, main=round(mean(b2), 8), xlab='', col = 'white', breaks = 40)
abline(v = mean(b2), col = 'blue', lwd=2, lty=2)
abline(v = Y['b2', ]['X2.5.'], col = 'red', lwd=2, lty=2)
abline(v = Y['b2', ]['X97.5.'], col = 'red', lwd=2, lty=2)
plot(acf(b2, lag.max = 100, plot = FALSE)[1:100])
par(mfrow=c(1,1)) 
################################################################################
############### Análise numérica
b = matrix( c( b0, b1, b2 ), ncol = 3 )
mcmcchain_b = coda::as.mcmc( b )
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_b = coda::geweke.diag( mcmcchain_b )
CD_b
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_b = coda::effectiveSize( mcmcchain_b )
IF_b = nrow(b) / N_eff_b
IF_b
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_b = round( apply( b, MARGIN = 2, FUN = sd) / sqrt( N_eff_b ), 5 )
mc_error_b
################################################################################

# h
chain_h = parameters$h
h_hat = apply(chain_h, 2, mean)
lim_inf = apply(chain_h, 2, quantile, probs=c(.025, .975))[1,]
lim_sup = apply(chain_h, 2, quantile, probs=c(.025, .975))[2,]
df = data.frame(x=1:T, ver = h, est = h_hat, 
                lim_inf = lim_inf, lim_sup = lim_sup)
g = ggplot(df[1:250, ]) 
g = g + geom_ribbon(aes(x, ymin = lim_inf, ymax = lim_sup), fill="gray80")
g = g + geom_line(aes(x=x, y = ver), color ='red')
g = g + geom_line(aes(x=x, y=est))
#g = g + geom_line(aes(x=x, y = lim_inf), linetype = 'dashed')
#g = g + geom_line(aes(x=x, y = lim_sup), linetype = 'dashed')
g
############### Análise numérica
mcmcchain_h = coda::as.mcmc( chain_h ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_h = coda::geweke.diag( mcmcchain_h )
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_h = coda::effectiveSize( mcmcchain_h )
IF_h = nrow(chain_h) / N_eff_h
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_h = round( apply( chain_h, MARGIN = 2, FUN = sd) / sqrt( N_eff_h ), 5 )

# plots
par( mfrow = c(1,3) )
plot( CD_h$z )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_h )
abline(h = 1)
plot( mc_error_h )
par( mfrow = c(1,1) )
###############################################################################
###############################################################################
############################## Model Selection
# construindo theta_hat e theta_draws
b_hat = data[4:6, 1]


theta_hat = c( b_hat, h_hat )
theta_draws = chain_b
theta_draws = rbind( theta_draws, chain_h )
theta_draws = as.matrix( theta_draws )
############### dic deviance information criterion:
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# theta = b, h
# data = c( y0, y )
# theta_hat = ( b_hat, h_hat )
# theta_draws = burned_lag
log_lik = function(theta_t, data){
  # função checada (13/04/23)
  T = length( data ) - 1
  
  b0_t = theta_t[1]
  b1_t = theta_t[2]
  b2_t = theta_t[3]
  h_t = theta_t[4:(T+3)]
  
  log_l = dnorm(data[2:( T + 1 )], 
                mean = b0 + b1 * data[1:T] + b2 * exp( h_t ), 
                sd = exp( 0.5 * h_t ), 
                log = TRUE)
  
  return( sum( log_l ) )
  
}
log_lik_i = function(data, theta_draws){
  # função checada (13/04/23)
  x = apply(X = theta_draws, MARGIN = 2, FUN = log_lik, data)
  
  return( x )
}
dic = function(data, theta_draws, theta_hat){
  
  pD = 2 * ( log_lik( theta_hat, data ) - mean( log_lik_i(data, theta_draws) ) )      
  DIC = - 2 * log_lik( theta_hat, data ) + 2 * pD  
  
  return( DIC )
} 

# calculando DIC
dic_svmn = dic( data = c(y0, y), 
                theta_draws = theta_draws, 
                theta_hat )

############### loo
lik = function(data_i, draws, data_, data_past){
  #data_ = ( y_{1}, y_{2}, ..., y_{T} )
  #data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  k = which( data_ == as.numeric( data_i ) )
  log_l = NULL
  N = ncol( draws )
  
  for(col in 1:N){
    
    b0_draws = draws[1, col]
    b1_draws = draws[2, col]
    b2_draws = draws[3, col]
    h_draws  = draws[3 + k, col]
    
    log_l[col] = dnorm(data_i, mean = b0_draws + b1_draws * data_past[k] + b2_draws * exp( h_draws ), 
                       sd = exp( 0.5 * h_draws ) )
    
  }
  
  return( log_l )
}

r_eff = loo::relative_eff(lik,
                          chain_id = rep(1, ncol( theta_draws ) ),
                          data = matrix( y , ncol = 1), 
                          draws = theta_draws,
                          data_ = y,
                          data_past = c( y0, y[1:(T-1)] )
                          #cores = getOption('mc.cores', 3)
)

# or set r_eff = NA
loo::loo(lik, 
         #r_eff = NA,
         r_eff = r_eff, 
         data = as.matrix( y ), 
         draws = theta_draws,
         data_ = y,
         data_past = c( y0, y[1:(T-1)] ),
         cores = getOption('mc.cores', 3)
)

############### waic
log_lik = function(data_i, draws, data_, data_past){
  #data_ = ( y_{1}, y_{2}, ..., y_{T} )
  #data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  k = which( data_ == as.numeric( data_i ) )
  log_l = NULL
  N = ncol( draws )
  
  for(col in 1:N){
    
    b0_draws = draws[1, col]
    b1_draws = draws[2, col]
    b2_draws = draws[3, col]
    h_draws  = draws[3 + k, col]
    
    log_l[col] = dnorm(data_i, mean = b0_draws + b1_draws * data_past[k] + b2_draws * exp( h_draws ), 
                       sd = exp( 0.5 * h_draws ), log = TRUE )
    
  }
  
  return( log_l )
}

loo::waic(log_lik, 
          data = matrix( y, ncol = 1 ), 
          draws = theta_draws,
          data_ = y,
          data_past = c( y0, y[1:(T-1)] )
)

