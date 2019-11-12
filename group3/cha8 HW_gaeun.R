#8.3
#(a)
##data
library(runjags)
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd('C:/Users/MYPC/Desktop')
school<-read.csv(file='school.csv', header=TRUE)

##MCMC
###graph
modelString <- "
model{

###likelihood
for (i in 1:N){
y[i] ~ dnorm(mu_j[school[i]], invsigma2)
}

###hyperprior
mu ~ dnorm(mu0, 1/g0^2)
invtau2 ~ dgamma(a_t, b_t)
tau <- sqrt(pow(invtau2, -1))

###prior
for (j in 1:J){
mu_j[j] ~ dnorm(mu, invtau2)
}
invsigma2 ~ dgamma(a_g, b_g)
sigma <- sqrt(pow(invsigma2, -1))
}
"

##Gibbs Sampler
y <- school %>% pull(hours)
school <- school %>% pull(schools)
N <- length(y)
J <- length(unique(school))

initsfunction <- function(chain){
  .RNG.seed <- c(1,2)[chain]
  .RNG.name <- c('base::Super-Duper',
                 'base::Wichmann-Hill')[chain]
  return(list(.RNG.seed=.RNG.seed,
              .RNG.name=.RNG.name))
}

the_data <- list('y'=y, 'school'=school, 'N'=N, 'J'=J,
                 'mu0'=7, 'g0'=5,
                 'a_t'=10, 'b_t'=2,
                 'a_g'=15, 'b_g'=2)


##JAGS
library(runjags)
runjags.options(silent.runjags=TRUE, silent.jags=TRUE)
posterior <- run.jags(modelString,
                      n.chains = 1,
                      data = the_data,
                      monitor = c('mu', 'tau', 'mu_j', 'sigma'),
                      adapt = 1000,
                      burnin = 5000,
                      sample = 1000,
                      thin = 1,
                      inits = initsfunction)
options(digits=2)
summary(posterior)

##anoter MCMC
library(dplyr)
library(tidyr)
library(ggplot2)

##data
schools.list=lapply(1:8,function(i){
  s.tbl=paste0('http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/school',i,'.dat')%>%
    url%>%read.table
  data.frame(school=i, hours=s.tbl[,1]%>%as.numeric)
})
schools.data<-do.call(rbind, schools.list)
Y<-schools.data

##prior
mu0<-7; g20<-5
t20<-10; eta0<-2
s20<-15; nu0<-2

##starting values
m<-length(unique(Y[,1]))
n<-sv<-ybar<-rep(NA,m)
for (j in 1:m)
{
  Yj<-Y[Y[,1]==j,2]
  ybar[j]<-mean(Yj)
  sv[j]<-var(Yj)
  n[j]<-length(Yj)
}
theta<-ybar; sigma2<-mean(sv)
mu<-mean(theta); tau2<-var(theta)

##setup MCMC
set.seed(123)
S<-1000
THETA<-matrix(nrow=S, ncol=m)
SMT<-matrix(nrow=S, ncol=3)

##MCMC algorithm
for(s in 1:S)
{
  ###sample new thetas
  for(j in 1:m)
  {
    vtheta<-1/(n[j]/sigma2+1/tau2)
    etheta<-vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j]<-rnorm(1,etheta,sqrt(vtheta))
  }
  ##sample new sigma2
  nun<-nu0+sum(n)
  ss<-nu0*s20
  for(j in 1:m){ss<-ss+sum((Y[Y[,1]==j,2]-theta[j])^2)}
  sigma2<-1/rgamma(1,nun/2,ss/2)
  ##sample new mu
  vmu<-1/(m/tau2+1/g20)
  emu<-vmu*(m*mean(theta)/tau2+mu0/g20)
  mu<-rnorm(1,emu,sqrt(vmu))
  ##sample a new tau2
  etam<-eta0+m
  ss<-eta0*t20+sum((theta-mu)^2)
  tau2<-1/rgamma(1,etam/2,ss/2)
  ##store results
  THETA[s,]<-theta
  SMT[s,]<-c(sigma2,mu,tau2)
}

##effective sample size
library(coda)
effectiveSize(SMT[,1])
effectiveSize(SMT[,2])
effectiveSize(SMT[,3])

#(b)
##95% C.I
quantile(SMT[,1], probs=c(0.025, 0.975))
quantile(SMT[,2], probs=c(0.025, 0.975))
quantile(SMT[,3], probs=c(0.025, 0.975))

##densities
install.packages("MCMCpack")
library(MCMCpack)

sigma2_prior = data.frame(
  value = seq(10, 22, by = 0.1),
  density = 1/dgamma(seq(10, 22, by = 0.1), nu0 / 2, nu0 * s20 / 2),
  variable = 'sigma2'
)
tau2_prior = data.frame(
  value = seq(0, 40, by = 0.1),
  density = 1/dgamma(seq(0, 40, by = 0.1), eta0 / 2, eta0 * t20 / 2),
  variable = 'tau2'
)
mu_prior = data.frame(
  value = seq(0, 12, by = 0.1),
  density = dnorm(seq(0, 12, by = 0.1), mu0, sqrt(g20)),
  variable = 'mu'
)
priors = rbind(sigma2_prior, tau2_prior, mu_prior)
priors$dist = 'prior'

smt.df=data.frame(SMT)
colnames(smt.df)=c('sigma2', 'mu', 'tau2')
smt.df$s=1:S
cut_size=10
smt.df = smt.df %>%
  tbl_df %>%
  mutate(scut = cut(s, breaks = cut_size)) %>%
  gather('variable', 'value', sigma2:tau2)
smt.df$dist = 'posterior'
ggplot(priors, aes(x = value, y = density, color = dist)) +
  geom_line() +
  geom_density(data = smt.df, mapping = aes(x = value, y = ..density..)) +
  facet_wrap(~ variable, scales = 'free')


#(c)
tau_draws<-as.mcmc(posterior, vars='tau')
sigma_draws<-as.mcmc(posterior, vars='sigma')
R = tau_draws^2/(tau_draws^2 + sigma_draws^2)
df<-data.frame(R = R)
ggplot(df, aes(x=R)) +  geom_density() + labs(title='Density of R')
# R이 1에 가까울수록 그룹들간의 variability가 크다

#(d)
theta7_6 = THETA[,7]<THETA[,6]
mean(theta7_6)
theta7 = (THETA[,7]<THETA[,-7])%>%apply(MARGIN=1,FUN=all)
mean(theta7)

#(e) -> error
relationship=data.frame(sample_average=ybar, post_exp=colMeans(THETA), 
                        school=1:length(ybar))
ggplot(relationship, aes(x=sample_average, y=post_exp, label=school)) +
  geom_text() + geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=mean(schools.raw[,'hours']), lty=2) +
  annotate('text', x=10, y=7.9, label=paste0("Pooled sample mean", round(mean(schools.raw[, 'hours']), 2)))+
  geom_hline(yintercept=mean(SMT[,'mu']), color='red') +
  annotate('text', x=10, y=7.4, label=paste0("Posterior exp. mu", round(mean(SMT[,'mu']), 2)),color='red')
