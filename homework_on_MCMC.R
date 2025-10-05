# library(ggplot2)
# library(MASS)
# library(MCMCpack)   
set.seed(343428)
rm(list = ls())

################################             ###################################
################################################################################
################################ ESERCIZIO 1 ###################################
################################################################################
###############################              ###################################


################################################################################
################################# Punto 1 ######################################
################################################################################

# Parametri

tau2 = 0.5
sigma2 = 5
phi = 3/10
beta = c(2,0.1)


# Creazione griglia da uniforme

n_points = 100 
coords = data.frame(long = rep(NA,n_points), lat = rep(NA,n_points))
for (i in 1:n_points)
{
  coords[i, ] <- c(runif(1, 0, 10), runif(1, 0, 10))
}


# Simulazione delle Y nei punti della griglia

m_s = beta[1] + beta[2]*coords[,1]
dist_mat = as.matrix(dist(coords))
cov_W = sigma2*(exp(-phi * dist_mat))
w_sample = t(chol(cov_W)) %*% matrix(rnorm(n_points), ncol = 1) + m_s  # simulazione di W
y_sample = w_sample + rnorm(n_points,0,sqrt(tau2))

################################################################################
################################# Punto 2 ######################################
################################################################################

# Creazione data frame utilizzato per il plot

data_sim = data.frame(long = coords$long, lat = coords$lat, w = w_sample, y = y_sample)
quantili <- quantile(data_sim$y, probs = c(0.25, 0.50, 0.75))

# Aggiunta della colonna "gruppo" utilizzata per assegnare colori diversi ai punti
data_sim$gruppo <- NA 

#Partizione delle osservazioni nei vari gruppi in base alla loro posizione rispetto ai quantili
for (i in 1:n_points) {
  y_i = data_sim$y[i]
  
  if (y_i <= quantili[1]) {
    data_sim$gruppo[i] <- "0-25% quantile"  
  } 
  else if (y_i <= quantili[2]) {
    data_sim$gruppo[i] <- "25-50% quantile"
  } 
  else if (y_i <= quantili[3]) {
    data_sim$gruppo[i] <- "50%-75% quantile"  
  } 
  else {
    data_sim$gruppo[i] <- "75%-100% quantile"
  }
}


# Scatterplot

dev.new()
ggplot(data_sim, aes(x = long, y = lat, color = gruppo)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("blue", "green", "orange", "red")) +  
  labs(x = "Longitudine", y = "Latitudine", color = "Gruppo") 

################################################################################
################################# Punto 3 ######################################
################################################################################

####################### Costruzioni dei dataset D_o e D_u ######################

n = 50  # numerosità di D_o

y_o_indici = sample(1:n_points, n)

# n osservazioni in y_o 
y_o = y_sample[y_o_indici]
w_o = w_sample[y_o_indici]

# Rispettive coordinate in D_o
D_o = data.frame(long = rep(NA,n), lat = rep(NA,n))
for (i in 1:n)
{
  D_o[i, ] <- c(coords$long[y_o_indici[i]], coords$lat[y_o_indici[i]])
}

# Altre osservazioni in y_u
y_u_indici = c()
for (i in 1:n_points)
{
  if (!(i %in% y_o_indici))
  {
    y_u_indici = c(y_u_indici,i)
  }
}

y_u = y_sample[y_u_indici]

# Rispettive coordinate in D_u
D_u = data.frame(long = rep(NA,n_points-n), lat = rep(NA,n_points-n))
for (i in 1:(n_points-n))
{
  D_u[i, ] <- c(coords$long[y_u_indici[i]], coords$lat[y_u_indici[i]])
}


####################### Distanza max e min nella griglia #######################

# Calcolo distanza massima e minima tra punti della griglia costruita al punto 1

min_distance <- Inf
max_distance <- -Inf

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    distance = sqrt((coords$long[j] - coords$long[i])^2 + (coords$lat[j] - coords$lat[i])^2)
    if (distance < min_distance) {
      min_distance = distance
    }
    if (distance > max_distance) {
      max_distance = distance
    }
  }
}


################## Costruzione matrice di covarianza di W ######################

# funzione che restituisce la matrice di covarianza di W utilizzata nel MCMC

matrix_cov <- function(X, phi, sigma2) {
  
  n = nrow(X)
  M = matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in i:n) {   
      dist = sqrt((X[i,1 ] - X[j,1])^2 + (X[i,2] - X[j,2])^2 )
      cov_value = sigma2 * exp(-phi * dist)
      M[i, j] = cov_value
      M[j, i] = cov_value  
    }
  }
  return(M)
}

####################### Funzione per algoritmo MCMC ############################

bayesian_mcmc_withBurnInAndThin <- function(W, Y, X, mean_prior_beta, m_beta, a_sigma , b_sigma , a_phi, b_phi, a_tau, b_tau, n_iter,
                                            adapt_batch = 50, adapt_a = 2000, adapt_b = 1000, adapt_alpha = 0.234, adapt_stop = 1000,
                                            burn_in_beta0, burn_in_beta1, burn_in_sigma2, burn_in_tau2, burn_in_phi,
                                            thin_beta0, thin_beta1, thin_sigma2, thin_tau2, thin_phi)
{
  
  # Numero di osservazioni
  n <- length(Y)
  
  
  # Inizializzazione vettori che contengono i campioni ad ogni iterazione
  
  beta_samples <- matrix(NA, n_iter, ncol(X))
  sigma2_samples <- numeric(n_iter)
  gamma_samples <- numeric(n_iter)
  tau2_samples <- numeric(n_iter)
  phi_samples <- numeric(n_iter)
  W_samples <- matrix(NA, n_iter, n)
  
  S = cbind(rep(1, n), X[,1])  # contiene una prima colonna di 1 e una seconda con la prima coordinata di X
  
  
  # Valori iniziali
  
  beta <- rep(0, ncol(X))  
  sigma2 <- 1
  tau2 <- 1
  phi <- runif(1,a_phi, b_phi)
  mat_cov_W <- matrix_cov(X, phi, sigma2)
  inv_mat_cov_W <- solve(mat_cov_W)
  m_s = S %*% beta
  
  # Parametri per la proposta adattiva:
  
  sd_proposal = 0.1 # Deviazione standard iniziale per la proposal
  somma_alpha = 0   # Contatore che tiene conto dei rapporti metropolis
  
  # Algoritmo MCMC:
  
  for (i in 1:n_iter) 
  {
    print("Iterazione numero:")
    print(i)
    
    # Aggiornamento beta 
    V_beta = solve(t(S) %*% inv_mat_cov_W %*% S + solve(m_beta))
    M_beta = V_beta %*% t(S) %*% inv_mat_cov_W %*% W
    beta = mvrnorm(1, M_beta, V_beta)     # Campionamento da una normale multivariata
    # Aggiornamento altre quantità relative a beta
    m_s = S %*% beta
    
    # Aggiornamento sigma
    C = sigma2 * inv_mat_cov_W
    qnt = t(W-m_s) %*% C %*% (W-m_s)
    shape_p = a_sigma + n/2
    rate_p = b_sigma + (1/2)*qnt   
    sigma2 = 1/rgamma(1,shape = shape_p, rate = rate_p)  # Campionamento da una inverse gamma
    # Aggiornamento altre quantità relative a sigma2
    mat_cov_W <- matrix_cov(X, phi, sigma2)
    inv_mat_cov_W <- solve(mat_cov_W)
    
    # Aggiornamento tau
    qnt = t(Y - W) %*% (Y - W)
    shape_post = a_tau + n/2
    rate_post = b_tau + (1/2)*qnt
    tau2 = 1/rgamma(1,shape = shape_post, rate = rate_post)  # Campionamento da una inverse gamma
    
    # Aggiornamento W
    V_w = solve(diag(nrow(inv_mat_cov_W))/tau2+ inv_mat_cov_W)
    M_w = V_w %*% (Y/tau2 + inv_mat_cov_W %*% m_s)
    W <- mvrnorm(1, M_w, V_w)                      # Campionamento da una normale multivariata
    W <- matrix(W, ncol = 1)
    
    # Aggiornamento phi tramite passo Metropolis
    eta = log((phi - a_phi)/(b_phi - phi))   # trasformazione in R
    eta_proposto <- rnorm(1, eta, sd_proposal)
    phi_proposto = (a_phi + b_phi*exp(eta_proposto))/(1 + exp(eta_proposto))  # trasformazione inversa
    
    # calcolo rapporto metropolis
    log_MH_ratio <- 0
    determ = determinant(mat_cov_W, logarithm = TRUE)$modulus
    mat_cov_W_prop <- matrix_cov(X, phi_proposto, sigma2)
    determ_prop = determinant(mat_cov_W_prop, logarithm = TRUE)$modulus
    inv_mat_cov_W_prop <- solve(mat_cov_W_prop)
    log_likeliwood_proposta = -(1/2)*(determ_prop + t(W-m_s) %*% inv_mat_cov_W_prop %*% (W-m_s))
    log_likeliwood_value_prec = -(1/2)*(determ + t(W-m_s) %*% inv_mat_cov_W %*% (W-m_s))
    log_prior_proposta = (eta_proposto - 2 * log(1 + exp(eta_proposto)))
    log_prior_value_prec = (eta - 2 * log(1 + exp(eta)))
    log_MH_ratio = log_likeliwood_proposta + log_prior_proposta - log_likeliwood_value_prec - log_prior_value_prec
    
    # Eventuale aggiornamento phi
    alpha = min(1, exp(log_MH_ratio)) 
    somma_alpha =  somma_alpha + alpha
    
    if (runif(1) < alpha)
    {
      # Accettazione della proposta
      phi = phi_proposto
      # Aggiornamento altre quantità relative a phi
      mat_cov_W <- matrix_cov(X, phi, sigma2)
      inv_mat_cov_W <- solve(mat_cov_W)
    }
    
    # Adattamento della varianza "sd_proposal"^2
    if (i %% adapt_batch == 0)
    {
      if (i < adapt_stop)
      {
        num_acc_phi_medio = somma_alpha / adapt_batch  # tasso medio di accettazione 
        sd_proposal = exp(log(sd_proposal) + adapt_a/(adapt_b + i) * (num_acc_phi_medio - adapt_alpha))
        somma_alpha = 0
      }
    }
    
    # Salvataggio campioni
    beta_samples[i, ] <- beta
    sigma2_samples[i] <- sigma2
    tau2_samples[i] <- tau2
    phi_samples[i] <- phi
    W_samples[i, ] <- W
  }
  
  # burn-in
  beta0_samplesB = beta_samples[burn_in_beta0:n_iter,1]
  beta1_samplesB = beta_samples[burn_in_beta1:n_iter,2]
  sigma2_samplesB = sigma2_samples[burn_in_sigma2:n_iter]
  tau2_samplesB = tau2_samples[burn_in_tau2:n_iter]
  phi_samplesB = phi_samples[burn_in_phi:n_iter]
  
  # thinning
  beta0_samplesBT = beta0_samplesB[seq(1, length(beta0_samplesB), by = thin_beta0)]
  beta1_samplesBT = beta1_samplesB[seq(1, length(beta1_samplesB), by = thin_beta1)]
  sigma2_samplesBT = sigma2_samplesB[seq(1, length(sigma2_samplesB), by = thin_sigma2)]
  tau2_samplesBT = tau2_samplesB[seq(1, length(tau2_samplesB), by = thin_tau2)]
  phi_samplesBT = phi_samplesB[seq(1, length(phi_samplesB), by = thin_phi)]
  
  # Risultati
  list(
    beta0_samples = beta0_samplesBT,
    beta1_samples = beta1_samplesBT,
    sigma2_samples =  sigma2_samplesBT,
    tau2_samples = tau2_samplesBT,
    phi_samples = phi_samplesBT,
    W_samples = W_samples
  )
}

########################### Definizione delle prior ############################

# Prior di beta_i: normale di media 0 e varianza 10000
mean_prior_beta = 0
matrix_cov_beta = matrix(c(10000, 0, 0, 10000), nrow = 2, ncol = 2)

# Prior di sigma_2, IG di parametri:
a_sigma = 1
b_sigma = 1

# Prior di tau_2, IG di parametri:
a_tau = 1
b_tau = 1

# Prior di phi, uniforme di parametri:
a_phi = 3/max_distance
b_phi = 3/min_distance

# Burn in e thinning 
burn_in_phi = 500 
burn_in_beta0 = 1
burn_in_beta1 = 1
burn_in_sigma2 = 500 
burn_in_tau2 = 500  
thin_beta0 = 1 
thin_beta1 = 1 
thin_sigma2 = 10  
thin_tau2 = 10   
thin_phi = 10   

n_iter = 10500     # numero di iterazioni del MCMC
x_o = unname(as.matrix(D_o))

post_samples <- bayesian_mcmc_withBurnInAndThin(W = w_o,Y = y_o, X = x_o, mean_prior_beta, matrix_cov_beta, a_sigma, b_sigma , a_phi, b_phi,  a_tau , b_tau , n_iter,  
                              adapt_batch = 50, adapt_a = 500, adapt_b = 1000, adapt_alpha = 0.3, adapt_stop = 1000,
                              burn_in_beta0, burn_in_beta1, burn_in_sigma2, burn_in_tau2, burn_in_phi,
                              thin_beta0, thin_beta1, thin_sigma2, thin_tau2, thin_phi)

# Risultato

dev.new()  
stimatore_media_beta_0 = mean(post_samples$beta0_samples)
plot(seq(1, length(post_samples$beta0_samples)), post_samples$beta0_samples, type="l", main="Campioni a posteriori di Beta0", lwd = 0.75)
abline(h = beta[1], col = "orange", lwd = 4, lty = 1)
abline(h = stimatore_media_beta_0, col = "green", lwd = 4, lty = 2)

dev.new()  
stimatore_media_beta_1 = mean(post_samples$beta1_samples)
plot(seq(1, length(post_samples$beta1_samples)), post_samples$beta1_samples, type="l", main="Campioni a posteriori di Beta1", lwd = 0.75)
abline(h = beta[2], col = "orange", lwd = 4, lty = 1)
abline(h = stimatore_media_beta_1, col = "green", lwd = 4, lty = 2)

dev.new() 
stimatore_media_sigma2 = mean(post_samples$sigma2_samples)
plot(seq(1, length(post_samples$sigma2_samples)), post_samples$sigma2_samples, type="l", main="Campioni a posteriori di Sigma2", lwd = 1)
abline(h = sigma2, col = "orange", lwd = 4, lty = 1)
abline(h = stimatore_media_sigma2, col = "green", lwd = 4, lty = 2)

dev.new()  
stimatore_media_tau2 = mean(post_samples$tau2_samples)
plot(seq(1, length(post_samples$tau2_samples)), post_samples$tau2_samples, type="l", main="Campioni a posteriori di Tau2", lwd = 1)
abline(h = tau2, col = "orange", lwd = 4, lty = 1)
abline(h = stimatore_media_tau2, col = "green", lwd = 4, lty = 2)

dev.new() 
stimatore_media_phi = mean(post_samples$phi_samples)
plot(seq(1, length(post_samples$phi_samples)), post_samples$phi_samples, type="l", main="Campioni a posteriori di Phi", lwd = 1)
abline(h = phi, col = "orange", lwd = 4, lty = 1)
abline(h = stimatore_media_phi, col = "green", lwd = 4, lty = 2)



# Studio thinning 
dev.new()
par(mfrow = c(2,1))
# Calcolo e grafico dell'ACF per beta0 e beta1
acf(post_samples$beta0_samples, lag.max = 30,  main = "Autocorrelazione beta0") 
acf(post_samples$beta1_samples,  lag.max = 30, main = "Autocorrelazione beta1") 
dev.new()
par(mfrow = c(3, 1))
# Calcolo e grafico dell'ACF per sigma2
acf(post_samples$sigma2_samples, lag.max = 30, main = "Autocorrelazione Sigma^2", cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
# Calcolo e grafico dell'ACF per tau2
acf(post_samples$tau2_samples, lag.max = 30, main = "Autocorrelazione Tau^2",cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
# Calcolo e grafico dell'ACF per phi
acf(post_samples$phi_samples, lag.max = 30, main = "Autocorrelazione Phi", , cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)


################################################################################
################################# Punto 4 ######################################
################################################################################

# Campioni a posteriori ottenuti nel punto 3
beta0_post = post_samples$beta0_samples
beta1_post = post_samples$beta1_samples
sigma2_post = post_samples$sigma2_samples
tau2_post = post_samples$tau2_samples
phi_post = post_samples$phi_samples

# Numero di campioni di epsilon(s) che verranno calcolati in ogni punto
N <- min(length(beta0_post), length(beta1_post), length(tau2_post))

# "M_epsilon" -> matrice che conterrà i campionia posteriori di espilon(s) con s in D_o:
# ciascuna colonna corrisponde ai campioni generati da epsilon(s) con s punto fissato in D_o 
# [M]_ij = campione i_esimo di epsilon(s_j)
M_epsilon <- matrix(0, N, nrow(D_o))

for (col in 1:nrow(D_o))
{
  for (row in 1:N) 
  {
    Y_s <- y_o[col]
    
    # Calcolo di un campione a posteriori da epsilon(s) 
    M_epsilon[row,col] <- Y_s - (beta0_post[row] + beta1_post[row] * D_o$long[col])
  }
}

# Media campionaria dei campioni a posteriori per ogni punto in D_o
media_samplesPost_eps <- apply(M_epsilon, 2, mean)

# Rappresentazione grafica rispetto ai quantili
data_sim2 = data.frame(long = D_o$long, lat = D_o$lat, eps = media_samplesPost_eps)
quantili <- quantile(data_sim2$eps, probs = c(0.25, 0.50, 0.75))
data_sim2$gruppo <- NA  
for (i in 1:nrow(data_sim2)) {
  eps_i = data_sim2$eps[i]
  
  if (eps_i <= quantili[1]) {
    data_sim2$gruppo[i] <- "Gruppo 0-25% quantile"  
  } 
  else if (eps_i <= quantili[2]) {
    data_sim2$gruppo[i] <- "Gruppo 25-50% quantile"
  } 
  else if (eps_i <= quantili[3]) {
    data_sim2$gruppo[i] <- "Gruppo 50-75% quantile"  
  } 
  else {
    data_sim2$gruppo[i] <- "Gruppo 75-100% quantile"  
  }
}
dev.new()
ggplot(data_sim2, aes(x = long, y = lat, color = gruppo)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("blue", "green", "orange", "red")) +  
  labs( x = "Longitudine", y = "Latitudine", color = "Gruppo")



################################################################################
################################# Punto 5 ######################################
################################################################################

n_2 = 20
D_u_star <- data.frame(long = rep(NA, n_2), lat = rep(NA, n_2))

# Estrazione dei primi 20 punti in D_u in  e salvataggio in D_u_star
for (i in 1:n_2) 
{
  D_u_star[i, ] <- c(coords$long[y_u_indici[i]], coords$lat[y_u_indici[i]])
}

# Estrazione delle 20 osservazioni in D_u_star 
y_u_star = y_u[1:n_2]

# Numero di campioni di Y(s) calcolati per ogni punto 
N <- min(length(beta0_post), length(beta1_post), length(phi_post), length(sigma2_post), length(tau2_post))

# "M_y_star" struttura dove salvo i campioni a posteriori di y(s) con s in D_u_star.
# "M_y_star" -> matrice che conterrà i campioni di y(s) con s appartenente a D_u_star:
# ciascuna colonna corrisponde ai campioni generati da y(s) con s punto fissato in D_o 
# [M]_ij = campione i_esimo di epsilon(s_j)
M_y_star <- matrix(0, N, nrow(D_u_star))

x = unname(as.matrix(coords))   # coords contiene tutti i punti della griglia
S = cbind(rep(1, n_points), x[,1])

for (row in 1:N)
{
  # Simulazione di W a posteriori (normale multivariata sui punti della griglia)
  cov_W_post = sigma2_post[row] *(exp(-phi_post[row] * dist_mat)) # Matrice di covarianza a posteriori
  media_W = beta0_post[1] + beta1_post[2]*coords[,1]   # Media di W a posteriori
  W_s = mvrnorm(1,media_W, cov_W_post)  # campione di W multivariato a posteriori
  
  for (col in 1:nrow(D_u_star))
  {
    # Campione della marginale estratto dal campione della congiunta
    W_s_i = W_s[y_u_indici[col]]
    # Simulazione di Y_s
    Y_s_i <- rnorm(1, W_s_i, sqrt(tau2_post[row]))  
    M_y_star[row, col] <- Y_s_i
  }
}

# Vettori dove sono salvati gli estremi dell'intervallo di credibilità di Y(s_i) con s_i fissato in D_u_star
lq_vec = rep(NA,n_2)
uq_vec = rep(NA,n_2)


# estraggo informazioni sulla distribuzione a posteriori di Y(s) tramite la funzione summary()
for (i in 1:n_2)
{
  print("Punto:")
  print(i)
  print(summary(M_y_star[,i]))
  cat("valore vero assunto da y:", y_u_star[i])
  print("lower bound intervallo di credibilità:")
  lq = quantile(M_y_star[,i], probs = 0.25)
  lq_vec[i] = lq["25%"] 
  print(lq_vec[i])
  print("upper bound intervallo di credibilità:")
  uq = quantile(M_y_star[,i], probs = 0.975)
  uq_vec[i] = uq["97.5%"]  
  print(uq_vec[i])
  print(" ")
}

data_sim2 <- D_u_star
data_sim2$y_true <- y_u_star

data_sim2$gruppo <- NA  
for (i in 1:nrow(D_u_star)) 
{
  y_i = y_u_star[i]
  if (y_i <= uq_vec[i] && y_i >= lq_vec[i] )
  {
    data_sim2$gruppo[i] <- "Compreso nell'intervallo di confidenza" 
  } 
  else 
  {
    data_sim2$gruppo[i] <- "Non compreso nell'intervallo di confidenza"  
  }
}

dev.new()
ggplot(data_sim2, aes(x = long, y = lat, color = gruppo)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("black", "red")) +  
  labs(x = "Longitudine", y = "Latitudine", color = "Gruppo") 



################################             ###################################
################################################################################
################################ ESERCIZIO 2 ###################################
################################################################################
###############################              ###################################


rgumbel <- function(n,a,b)
{
  return(-log(-log(runif(n,a,b))))
}

# Parametri
n <- 200
pi_k <- c(1/3, 1/3, 1/3)  # Probabilità P(Z_i = k)
lambda_k <- c(1, 10, 25)  # Parametri di Poisson per ciascun gruppo

# Campionamento diretto di Y
z_sample <- sample(1:3, n, replace = TRUE, prob = pi_k)
y_sample <- rpois(n, lambda_k[z_sample])

# Rappresentazione grafica: istogrammi di frequenza per ciascun gruppo
group1 <- y_sample[which(z_sample== 1)]
group2 <- y_sample[which(z_sample == 2)]
group3 <- y_sample[which(z_sample== 3)]
breaks <- seq(0, max(y_sample) + 1, by = 1)
hist1 <- hist(group1,breaks = breaks, plot = FALSE)
hist2 <- hist(group2,breaks = breaks, plot = FALSE)
hist3 <- hist(group3, breaks = breaks,plot = FALSE)
max_freq <- max(c(hist1$counts, hist2$counts, hist3$counts))

dev.new()
hist(group1, breaks = breaks, main="Istogrammi di frequenza di Y_i|Z_i = k",
     xlab = "Valori", ylab = "Frequenza", xlim = c(0, max(y_sample)+1), 
     ylim = c(0, max_freq), col= "lightblue", border="black")
hist(group2, breaks= breaks, col = "lightgreen", border="black",add= TRUE)
hist(group3, breaks = breaks, col ="lightpink", border="black", add = TRUE)
legend("topright", legend = c("λ = 1","λ = 10", "λ = 25"),
       fill=c("lightblue", "lightgreen", "lightpink"), border="black")


# MCMC

bayesian_mcmc <- function(Z,Y,a_lambda, b_lambda, vec_alpha, n_iter)
{
  # Numero di osservazioni
  n <- length(Y)
  
  # Inizializzazione vettori che conterranno i campioni ad ogni iterazione
  lambda_samples <- matrix(NA, n_iter, 3)
  Z_samples <- matrix(NA, n_iter, n)
  pi_samples <- matrix(NA, n_iter, 3)
  
  # Valori iniziali
  lambda_vec = c(1,1,1) # vettore con i lambda_k
  pi_vec = rdirichlet(1, c(1,1,1))  # vettore con pi_k
  I_1 = which(Z == 1) # indici di Z tc Z == 1
  I_2 = which(Z == 2)
  I_3 = which(Z == 3)
  n_vec = c(length(I_1),length(I_2), length(I_3))
  
  for (iter in 1:n_iter)
  {
    print(iter)
    
    # Aggiornamento lambda_1
    shape_full = a_lambda + sum(Y[I_1])
    rate_full = b_lambda + n_vec[1]
    lambda_vec[1] = rgamma(1, shape = shape_full  , rate = rate_full)
    
    # Aggiornamento lambda_2
    shape_full = a_lambda + sum(Y[I_2])
    rate_full = b_lambda + n_vec[2]
    lambda_vec[2] = rgamma(1, shape = shape_full  , rate = rate_full)
    
    # Aggiornamento lambda_3
    shape_full = a_lambda + sum(Y[I_3])
    rate_full = b_lambda + n_vec[3]
    lambda_vec[3] = rgamma(1, shape = shape_full  , rate = rate_full)
    
    # Aggiornamento Z
    for (i in 1:length(Z))
    {
      # Aggiornamento di z_i (tramite Gumbel trick)
      y_i = y_sample[i]
      g_vec = dpois(y_i,lambda_vec,log=TRUE)
      gumb_samples = rgumbel(3,0,1) 
      qnt = log(pi_vec) +  g_vec + gumb_samples
      z_i = which.max(qnt);
      Z[i] = z_i
    }
    
    # Aggiornamento quantità relative a Z
    I_1 = which(Z == 1)
    I_2 = which(Z == 2)
    I_3 = which(Z == 3)
    n_vec = c(length(I_1),length(I_2), length(I_3))
    
    # Aggiornamento pi_vec
    vec_alpha_full_cond = vec_alpha + n_vec
    pi_vec = rdirichlet(1, vec_alpha_full_cond)
    pi_vec = pi_vec[1, ]
    
    # Salvataggio campioni
    lambda_samples[iter, ] <- lambda_vec
    Z_samples[iter, ] <- Z
    pi_samples[iter, ] <- pi_vec
  }
  
  # Risultati
  list(
    lambda1_samples = lambda_samples[1: n_iter,1],
    lambda2_samples = lambda_samples[1: n_iter,2],
    lambda3_samples = lambda_samples[1: n_iter,3],
    Z_samples = Z_samples,
    pi1_samples =  pi_samples[,1],
    pi2_samples =  pi_samples[,2],
    pi3_samples =  pi_samples[,3]
  )
}


# Prior di lambda i: Gamma di parametri (1,0.01)
a_gamma = 1 
b_gamma = 0.01 

# Prior della distribuzione pi: Dirichlet di parametri (1,1,1)
vec_alpha = c(1,1,1)

n_iter = 3000    # numero iterazioni MCMC

post_samples <- bayesian_mcmc(Z = z_sample,Y = y_sample, a_lambda = a_gamma, b_lambda = b_gamma, vec_alpha = vec_alpha, n_iter)


# Si vedono le catene

dev.new()
par(mfrow = c(3, 1))  

# Lambda_1 
plot(seq(1, length(post_samples$lambda1_samples)), post_samples$lambda1_samples, type = "l", main = "Lambda_1",
     xlab = "Iterations", ylab = "Campioni di Lambda_1")
abline(h = lambda_k[1], col = "orange", lwd = 3, lty = 1)
abline(h = mean(post_samples$lambda1_samples), 
       col = "green", lwd = 3, lty = 2)  

# Lambda_2 
plot(seq(1, length(post_samples$lambda2_samples)), post_samples$lambda2_samples, type = "l", main = "Lambda_2",
     xlab = "Iterations", ylab = "Campioni di Lambda_2")
abline(h = lambda_k[2], col = "orange", lwd = 3, lty = 1)  
abline(h = mean(post_samples$lambda2_samples), 
       col = "green", lwd = 3, lty = 2)  

# Lambda_3 Posterior
plot(seq(1, length(post_samples$lambda3_samples)), post_samples$lambda3_samples, type = "l", main = "Lambda_3",
     xlab = "Iterations", ylab = "Campioni di Lambda_3")
abline(h = lambda_k[3], col = "orange", lwd = 3, lty = 1)  
abline(h = mean(post_samples$lambda3_samples), 
       col = "green", lwd = 3, lty = 2)  

dev.new()
par(mfrow = c(3, 1))  

# pi_1
plot(seq(1, length(post_samples$pi1_samples)), post_samples$pi1_samples, type = "l", main = "pi_1",
     xlab = "Iterations", ylab = "Campioni di pi_1 ")
abline(h = pi_k[1], col = "orange", lwd = 3, lty = 1)  
abline(h = mean(post_samples$pi1_samples), 
       col = "green", lwd = 3, lty = 2)  

# pi_2
plot(seq(1, length(post_samples$pi2_samples)), post_samples$pi2_samples, type = "l", main = "pi_2",
     xlab = "Iterations", ylab = "Campioni di pi_2")
abline(h = pi_k[2], col = "orange", lwd = 3, lty = 1)  
abline(h = mean(post_samples$pi2_samples), 
       col = "green", lwd = 3, lty = 2)  

# pi_3
plot(seq(1, length(post_samples$pi3_samples)), post_samples$pi3_samples, type = "l", main = "pi_3",
     xlab = "Iterations", ylab = "Campioni di pi_3")
abline(h = pi_k[3], col = "orange", lwd = 3, lty = 1)  
abline(h = mean(post_samples$pi3_samples), 
       col = "green", lwd = 3, lty = 2)


