## Data
school1 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school1.dat')
school2 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school2.dat')
school3 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school3.dat')
school4 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school4.dat')
school5 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school5.dat')
school6 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school6.dat')
school7 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school7.dat')
school8 = scan('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school8.dat')

## Data transformation
library(ggplot2)
library(tidyr)
library(dplyr)

schools.list = lapply(1:8, function(i) {
  s.tbl = paste0('http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/school', i, '.dat') %>%
    url %>%
    read.table
  
  data.frame(
    school = i,
    hours = s.tbl[, 1] %>% as.numeric
  )
})

schools.raw = do.call(rbind, schools.list)
y = schools.raw


## rjag..로 try
library(rjags)
library(runjags)

school.list <- list(school1,school2,school3,school4,school5,school6,school7, school8)
n = as.numeric(lapply(school.list, length))
ybar = as.numeric(lapply(school.list, mean))
var <- as.numeric(lapply(school.list, var))
y <- school[,1]

theta = ybar
s2 <- mean(var)
mu <- mean(theta)
tau2 <- var(theta)

## Modeling
modelString <- "model{
  ## likelihood 
  for (i in 1:N){
    y[i] ~ dnorm(mu_j[school[i]], invsgima2)
  }
  
  
  ## hyperprior
  mu ~ dnorm(mu0, 1/g20)
  invtau2 ~ dgamma(a_t, b_t)
  tau <- sqrt(pow(invtau2, -1))
  
  ## prior
  for(j in 1:J){
  mu_j[j] ~ dnorm(mu, invtau2)
  }
  invsigma2 ~ dgamma(a_g,b_g)
  sigma <- sqrt(pow(invsigma2, -1))
  }
"

y <- school.data%>% pull(hours)
school <- school.data %>% pull(school)
N <- length(y)
J <- length(unique(school))

initsfunction <- function(chain){
  .RNG.seed <- c(1,2)[chain]
  .RNG.name <- c('base::Super-Duper',
                 'base::Wichmann-Hill')[chain]
  return(list(.RNG.seed=.RNG.seed,
              .RNG.name=.RNG.name))
}

the_data <- list('y'=y, 'school'=school, 'N'=N,'J'=J,
                 'mu0' = 7, 'g20'= 5, 'a_t'=10, 'b_t'=2, 'a_g'= 15, 'b_g'=2 )


## Prior
mu0 <- 7
g20 <- 5
t20 <- 10
eta0 <- 2
s20 <- 15
nu0<- 2

m<-length(unique(y[,1]))
n<-sv<-ybar<-rep(NA,m)

for (j in 1:m) {
  y_j<-y[y[,1]==j,2]
  ybar[j]<-mean(y_j)
  sv[j]<-var(y_j)
  n[j]<-length(y_j)
}

theta<-ybar
sigma2<-mean(sv)
mu<-mean(theta)
tau2<-var(theta)

## MCMC
set.seed(730)
S<-1000
Theta<-matrix(nrow=S, ncol=m)
SMT<-matrix(nrow=S, ncol=3)


for(s in 1:S)
{
  # sample thetas
  for(j in 1:m){
    vtheta <- 1/(n[j]/sigma2+1/tau2)
    etheta <- vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j] <- rnorm(1,etheta,sqrt(vtheta))
  }
  # sample sigma2
  nun<-nu0+sum(n)
  ss<-nu0*s20
  for(j in 1:m){
    ss <- ss+sum((y[y[,1]==j,2]-theta[j])^2)
    }
  sigma2 <-1/rgamma(1,nun/2,ss/2)
  
  # sample mu
  vmu <-1/(m/tau2+1/g20)
  emu <-vmu*(m*mean(theta)/tau2+mu0/g20)
  mu <- rnorm(1,emu,sqrt(vmu))
  
  # sample tau2
  etam <- eta0+m
  ss <- eta0*t20+sum((theta-mu)^2)
  tau2 <- 1/rgamma(1,etam/2,ss/2)
  
  # result
  Theta[s,] <-theta
  SMT[s,] <-c(sigma2,mu,tau2)
}

### (a)
library(coda)
effectiveSize(SMT[,1])
effectiveSize(SMT[,2])
effectiveSize(SMT[,3])

### (b)

# 95% Confidence Regions of sigma2, mu, tau2
t(apply(SMT, 2, quantile, probs=c(0.025,0.975)))

# compare the posterior to the prior
library(MCMCpack)

sigma2_prior = data.frame(
  value = seq(10, 22, by = 0.1),
  density = 1/dgamma(seq(10, 22, by=0.1), nu0/2, nu0*s20/2),
  variable = 'sigma2'
)
tau2_prior = data.frame(
  value = seq(0, 40, by = 0.1),
  density = 1/dgamma(seq(0, 40, by=0.1), eta0/2, eta0*t20/2),
  variable = 'tau2'
)
mu_prior = data.frame(
  value = seq(0, 12, by = 0.1),
  density = dnorm(seq(0, 12, by=0.1), mu0, sqrt(g20)),
  variable = 'mu'
)

priors = rbind(sigma2_prior, tau2_prior, mu_prior)

smt.df=data.frame(SMT)
colnames(smt.df)=c('sigma2','mu','tau2')
smt.df$s=1:S
cut_size=10
smt.df = smt.df %>%
  tbl_df %>%
  mutate(scut = cut(s, breaks = cut_size)) %>%
  gather('variable', 'value', sigma2:tau2)

priors$dist = 'prior'
smt.df$dist = 'posterior'

ggplot(priors, aes(x = value, y = density, color = dist)) +
  geom_line() +
  geom_density(data = smt.df, mapping = aes(x = value, y = ..density..)) +
  facet_wrap(~ variable, scales = 'free')


#(c)
tau_pos <-as.mcmc(posterior, vars='tau')
sigma_pos <-as.mcmc(posterior, vars='sigma')
R_prior = data.frame(tau_pos^2/(tau_pos^2 + sigma_pos^2))
R_posterior = data.frame(SMT[,'tau2']/SMT[,'tau2']+SMT[,'sigma2'])
ggplot(R_prior, aes(x=value, y=..density..)) +
  geom_density(data=R_prior)+
  geom_density(data=R_posterior)

#(d)
theta7smaller6 = Theta[,7]<Theta[,6]
mean(theta7smaller6)
theta7smallest = (Theta[,7]<Theta[,-7]) %>% apply(MARGIN=1,FUN=all)
mean(theta7smallest)

#(e) -> error
relationship=data.frame(
  sample_average=ybar, 
  post_exp=colMeans(Theta), 
  school=1:length(ybar)
  )

ggplot(relationship, aes(x=sample_average, y=post_exp, label=school)) +
  geom_text() + 
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=mean(schools.raw[,'hours']), lty=2) +
  annotate('text', x=10, y=7.9, label=paste0("Pooled sample mean", round(mean(schools.raw[, 'hours']), 2)))+
  geom_hline(yintercept=mean(SMT[,'mu']), color='red') +
  annotate('text', x=10, y=7.4, label=paste0("Posterior exp. mu", round(mean(SMT[,'mu']), 2)),color='red')
