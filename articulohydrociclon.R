############ Articulo Hidrociclon ###########

# Cargando paquetes
library(nlme)
library(Hmisc)
library(lattice)
library(rstan)
library(loo)
require(shinystan)
require(bayesplot)
library(readxl) # Leer excel
library(dplyr)

# Leer datos

datos <- read.csv("https://raw.githubusercontent.com/EmilioMonardes/Hydrocyclone/main/datos.csv", sep = ";")

##  Figura 

dados$ID=as.factor(dados$ID)
p1=ggplot(dados,aes(x=SIZE,y=Acum..Pasante,group=ID))+
    geom_line(aes(color=ID)) + xlab(
        expression(paste("Granularidad (",mu,"m)"))) + 
    ylab("Acumulado Pasante")+ theme(legend.position="none")

p2=ggplot(dados,aes(x=log(SIZE),y=Acum..Pasante,group=ID))+
    geom_line(aes(color=ID)) + xlab(
        expression(paste("Granularidad [log(",mu,"m)]"))) + 
    ylab("")

cowplot::plot_grid(p1,p2)



#Problemas en 

#m12 cambie la tolerancia
#m13 revisar algo del algoritmo
#m9 problemas
#m11

# Load packages                           
library(nlme)
library(Hmisc)
library(lattice)

# Load functions for elliptical models

source("https://raw.githubusercontent.com/EmilioMonardes/Hydrocyclone/main/Elliptical_functions.r")

# Choose elliptical model & nonlinear model
modelo <- 'normal'
funcao = 'm7'

# Number of observations by group/subject
m <- c(50,50,50,50,50,50,40,50,50,49,35,47,40,42,50)
mm <- m
ma = c(1,(cumsum(m)+1))

y<- dados$Y
x<- dados$X
dados$x <- x
n = length(y)
p<- 3  
ii = unique(dados$ID)

# Nonlinear function 
f<-function(t,x,funcao){
    beta1 = t[1];
    beta2 = t[2];
    beta3 = t[3];
    switch(funcao,
           'logistica' = {beta1/(1+exp(-(x-beta2)/beta3))},
           'Gompertz' = {
               .a = beta1 * exp(-beta2 * exp(-beta3 * x)); 
               .a[is.na(.a)] <- 0;
               .a;},
           'm2' = {100 / (1 + (beta1 / x)^beta2) ^ beta3},
           'm3' = {beta1 + ((100 - beta1)/(1 + (beta2 / x) ^ beta3) ^ ((1 - 1 / beta3)))},
           'm4' = {100 - ((100 - beta1) * exp((-(beta2*x)^beta3)) )},
           'm5' = {(beta3 * beta1 ^ beta2) / (beta1 ^ beta2 + x ^ beta2)},
           'm6' = {beta1 * (1 - ((((beta1 - beta2) / beta1) * exp(-beta3 * x))))},
           'm7' = {beta1 * (1 - exp(- beta2 * (x-beta3)))},
           'm8' = {beta1 * (1-(((beta1 - beta2) / beta1) * exp(-beta3 * x)))},
           'm9' = {100/(1+(100/ beta1 - 1) * exp(-beta2 * x^beta3))},
           'm10' = {100-((100-beta1)/(1+(beta2 * x)^beta3))},
           'm11' = {100 * (1 + (beta1 - 1) * exp(- beta2 *(x-beta3))) ^ (1/(1- beta1))},
           'm12' = {100/(1 + (log (200 /x) / log (200)/beta2 )) ^ beta1},
           'm16' = {beta3 *(x/ beta1)^beta2},
           'm14' = {100 * (1-exp(-(x / beta1) ^ beta2))}
           
    )
}

# Gradient matrix
J<-function(t,x,funcao){
    beta1 = t[1];
    beta2 = t[2];
    beta3 = t[3];
    switch(funcao,
           'logistica'=
               {
                   J1<- 1/(1+exp(-(x-beta2)/beta3));
                   J2<- -beta1*exp(-(x-beta2)/beta3)/((1+exp(-(x-beta2)/beta3))^2*beta3);
                   J3<- -beta1*(x-beta2)*exp(-(x-beta2)/beta3)/((1+exp(-(x-beta2)/beta3))^2*beta3^2);
                   cbind(J1,J2,J3);
               },
           'Gompertz' = 
               {
                   #f<-function(beta1,beta2,beta3,x){ beta1 * exp(-beta2 * exp(-beta3 * x))}
                   #Deriv(f,'beta1')
                   J1 = exp(-(beta2 * exp(-(beta3 * x))));
                   
                   #Deriv(f,'beta2')
                   .e1 <-  .e2 <- exp(-(beta3 * x));
                   J2 <- -(beta1 * exp(-(beta2 * .e2)) * .e2);
                   
                   #Deriv(f,'beta3')
                   .e2 <- exp(-(beta3 * x));
                   J3 <- beta1 * beta2 * x * exp(-(beta2 * .e2)) * .e2;    
                   cbind(J1,J2,J3);},
           'm2'= { #Havercamp
               #f function <- 100 / (1 + (beta1 / x)^beta2) ^ beta3; 
               J1<- -(100 * beta2 * beta3 * (beta1 / x)^beta2 * ((beta1 / x)^beta2 + 1)^(-beta3 - 1)) / beta1;
               J2<- -100 * (beta1 / x)^beta2 * log(beta1 / x) * ((beta1 / x)^beta2 + 1)^(-beta3 - 1) * beta3;
               J3<- -(100 * log((beta1 / x)^beta2 + 1)) / ((beta1 / x)^beta2 + 1)^beta3;
               cbind(J1,J2,J3);
           },
           'm3'= { #LIMA
               #f funcition <- beta1 + ((100 - beta1)/(1+ (beta2 / x)^beta3)^{(1 - 1/ beta3)}); 
               J1<- 1 - ((beta2 / x)^beta3 + 1)^(1 / beta3 - 1);
               J2<- ((100 - beta1) * (1 / beta3 - 1) * beta3 * (beta2 / x)^beta3 * ((beta2 / x)^beta3 + 1)^(1 / beta3 - 2)) / beta2;
               J3<- (100 - beta1) * ((beta2 / x)^beta3 + 1)^(1 / beta3 - 1) * (((beta2 / x)^beta3 * log(beta2 / x) * (1 / beta3 - 1)) / ((beta2 / x)^beta3 + 1) - log((beta2 / x)^beta3 + 1) / beta3^2);
               cbind(J1,J2,J3);
           },
           'm4'= { #Weibull
               #f function <- 100 - ((100 - beta1) * exp((-(beta2 * x)^beta3))); 
               J1<- exp(-(x * beta2)^beta3);
               J2<- ((100 - beta1) * beta3 * (x * beta2)^beta3 * exp(-(x * beta2)^beta3)) / beta2;
               J3<- (100 - beta1) * (x * beta2) ^ beta3 * log(x * beta2) * exp(-(x * beta2) ^ beta3);
               cbind(J1,J2,J3);
           },
           'm5'= { #Hill 
               #f funcition <- (beta3 * beta1 ^ beta2) / (beta1 ^ beta2 + x ^ beta2); 
               J1<- (x ^ beta2 * beta2 * beta3 * beta1 ^ (beta2 - 1)) / (beta1 ^ beta2 + x ^ beta2)^2;
               J2<- (x ^ beta2 * beta1 ^ beta2 * (log(beta1) - log(x)) * beta3) / (beta1 ^ beta2 + x ^ beta2)^2;
               J3<- beta1^ beta2 / (beta1 ^ beta2 + x ^ beta2);
               cbind(J1,J2,J3);
           },
           
           'm6'= { #Brody
               #f ; beta1 * (1 - ((((beta1 - beta2) / beta1) * exp(-beta3 * x)))
               J1<- exp(-x * beta3) * (exp(x * beta3) - 1);
               J2<- exp(-x * beta3);
               J3<- x * (beta1 - beta2) * exp(-x * beta3);
               cbind(J1,J2,J3);
           },
           
           'm7'= { #Bertalanffy 
               #f ; beta1 * (1 - exp(- beta2 * (x-beta3)))
               J1<- 1 - exp(-beta2 * (x - beta3));
               J2<- -beta1 * (beta3 - x) * exp(-(x - beta3) * beta2);
               J3<- -beta1 * beta2 * exp(-beta2 * (x - beta3));
               cbind(J1,J2,J3);
           },
           
           'm8'= { #Monomolecular 
               #f ; beta1 * (1-(((beta1 - beta2) / beta1) * exp(-beta3 * x)))
               J1<- exp(-x * beta3) * (exp(x * beta3) - 1);
               J2<- exp(-x * beta3);
               J3<- x * (beta1 - beta2) * exp(-x * beta3);
               cbind(J1,J2,J3);
           },
           
           'm9'= { #Skaggs
               #f ; 100/(1+(100/ beta1 - 1) * exp(-beta2 * x^beta3))
               J1<- (10000 * exp(x^beta3 * beta2)) / ((exp(x ^ beta3 * beta2) - 1) * beta1 + 100)^2;
               J2<- -(100 * x ^ beta3 * (beta1 - 100) * beta1 * exp(x ^ beta3 * beta2)) / (beta1 * exp(x ^ beta3 * beta2) - beta1 + 100)^2;
               J3<- -(100 * x ^ beta3 * log(x) * (beta1 - 100) * beta1 * beta2 * exp(x ^ beta3 * beta2)) / (beta1 * exp(x ^ beta3 * beta2) - beta1 + 100)^2;
               cbind(J1,J2,J3);
           },
           
           'm10'= { #Morganetal
               #f ; 100-((100-beta1)/(1+(beta2 * x)^beta3))
               J1<- 1 / ((x * beta2) ^ beta3 + 1);
               J2<- ((100 - beta1) * beta3 * (x * beta2) ^ beta3) / (beta2 * ((x * beta2)^beta3 + 1)^2);
               J3<- ((100 - beta1) * (x * beta2) ^ beta3 * log(x * beta2)) / ((x * beta2)^beta3 + 1)^2;
               cbind(J1,J2,J3);
           },
           
           'm11'= { #Richards
               #f ; 100 * (1 + (beta1 - 1) * exp(- beta2 *(x-beta3))) ^ (1/(1- beta1))
               J1<- 100 * (exp(-beta2 * (x - beta3)) / ((exp(-beta2 * (x - beta3)) * (beta1 - 1) + 1) * (1 - beta1)) + log(exp(-beta2 * (x - beta3)) * (beta1 - 1) + 1) / (1 - beta1)^2) * (exp(-beta2 * (x - beta3)) * (beta1 - 1) + 1)^(1 / (1 - beta1));
               J2<- -(100 * (beta3 - x) * ((beta1 - 1) * exp(-(x - beta3) * beta2) + 1)^(1 / (1 - beta1))) / (exp((x - beta3) * beta2) + beta1 - 1);
               J3<- -(100 * beta2 * ((beta1 - 1) * exp(-beta2 * (x - beta3)) + 1)^(1 / (1 - beta1))) / (exp(beta2 * (x - beta3)) + beta1 - 1);
               cbind(J1,J2,J3);
           },
           
           'm12'= { #Swebrec
               #f ; 100/(1 + (log (200 /x) / log (200)/beta2 )) ^ beta1
               J1<- -(log(log(200 / x) / (log(200) * beta2) + 1) * beta3) / (log(200 / x) / (log(200) * beta2) + 1)^beta1;
               J2<- (log(200 / x) * beta1 * beta3 * (log(200 / x) / (log(200) * beta2) + 1)^(-beta1 - 1)) / (log(200) * beta2^2);
               J3<- 1 / (log(200 / x) / (log(200) * beta2) + 1) ^ beta1;
               cbind(J1,J2,J3);
           },
           
           'm13'= { #Bass
               #f ; (beta1 * ((1-exp(-(beta2+beta3)* x))))/((1+ (beta3 / beta2) * exp(-(beta2+ beta3) * x)))
               J1<- (beta2 * (exp(x * (beta3 + beta2)) - 1)) / (beta2 * exp(x * (beta3 + beta2)) + beta3);
               J2<- (beta1 * ((x * beta2^2 + x * beta3 * beta2 + beta3) * exp(x * (beta2 + beta3)) - beta3)) / (beta2 * exp(x * (beta2 + beta3)) + beta3)^2;
               J3<- (beta1 * beta2 * ((x * beta3 + x * beta2 - 1) * exp(x * (beta3 + beta2)) + 1)) / (beta2 * exp(x * (beta3 + beta2)) + beta3)^2;
               cbind(J1,J2,J3);
           },
           'm14'= { #Rosin-Rammler
               #f ; beta3 * (1-exp(-(x / beta1) ^ beta2))
               J1<- -(beta2 * beta3 * (x / beta1)^beta2 * exp(-(x / beta1)^beta2)) / beta1;
               J2<- log(x / beta1) * (x / beta1)^beta2 * exp(-(x / beta1)^beta2) * beta3;
               J3<- 1 - exp(-(x / beta1)^beta2);
               cbind(J1,J2,J3);
           },
           'm16'= { #Gaudin-Schuhmann
               #f ; beta3 * (x / beta1) ^ beta2
               J1<- -(beta2 * beta3 * (x / beta1)^ beta2) / beta1;
               J2<- log(x / beta1) * beta3 * (x / beta1)^beta2;
               J3<- (x / beta1)^beta2;
               cbind(J1,J2,J3);
           },
           'm17'= { #Fredlund
               #f ; 100 / log(exp(1) + ((( beta1 / x) ^ {beta2})^{beta3})) * (1-(log(1+ min(x)/{x})/log(1+ min(x)/0.00001  ))^7)
               J1<- -(beta2 * beta3 * (x / beta1)^ beta2) / beta1;
               J2<- log(x / beta1) * (x / beta1)^beta2 * beta3;
               J3<- (x / beta1)^beta2;
               cbind(J1,J2,J3);
           }
    )
}

# Initial model
modelo_inicial <- 
    switch(funcao,
           'logistica'= {
               nls(y ~ K / (1 + exp(-b * (x - M))), 
                   start = list(K = max(y), b = 0.5, M = median(x)))},
           'Gompertz' = {nls(y ~ a * exp(-b * exp(-c * log(x))),
                             start = list(a = max(y), b = 10, c = 0.9),
                             control = nls.control(maxiter = 100, tol = 1e-6))},
           'm2' = {nls(y ~ 100 / (1 + (a / x) ^ b) ^ c,
                       start = list(a = 5, b = 2, c = 2),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm3' = {nls(y ~ a + ((100 - a)/(1+ (b / x) ^ c)^((1 - 1 / c))),
                       start = list(a = -1, b = 20, c = 2),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm4' = {nls(y ~ 100 - ((100 - a) * exp((-(b * x) ^ c))),
                       start = list(a = min(y), b = 0.029, c = 0.8),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm5' = {nls(y ~ (c * a ^ b) / (a ^ b + x ^ b),
                       start = list(a = 3, b = 2, c = max(y)),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm6' = {nls(y ~ a * (1 - ((((a - b) / a) * exp(-c * x)))),
                       start = list(a = max(y), b = min(y), c = 0.002),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm7' = {nls(y ~ a * (1 - exp(- b * (x - c))),
                       start = list(a = max(y), b = 0.01, c = min(x)),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm8' = {nls(y ~ a * (1-(((a - b) / a) * exp(-c * x))),
                       start = list(a = max(y), b = min(y), c = 0.001),
                       control = nls.control(maxiter = 100, tol = 1e-6))},
           'm9' = {nls(y ~ 100/(1+(100/ a - 1) * exp(-b * x^c)),
                       start = list(a = 2, b = 4, c = 0.1),
                       control = nls.control(maxiter = 300, tol = 1e-6))},
           'm10' = {nls(y ~ 100-((100-a)/(1+(b * x)^c)),
                        start = list(a = min(y), b = 0.1, c = 0.9),
                        control = nls.control(maxiter = 100, tol = 1e-6))},
           'm11' = {nls(y ~ 100 * (1 + (a - 1) * exp(- b *(x-c))) ^ (1/(1 - a)),
                        start = list(a = 0.2, b = 0.1, c = 10),
                        control = nls.control(maxiter = 100, tol = 1e-6))},
           'm12' = {nls(y ~ 100/(1 + (log (200 /x) / log (200)/b )) ^ a,
                        start = list(a = 10, b = 5),
                        control = nls.control(maxiter = 100, tol = 0.1))},
           'm13' = {nls(y ~ (beta1 * ((1-exp(-(beta2+beta3)* x))))/((1+ (beta3 / beta2) * exp(-(beta2+ beta3) * x))),
                        start = list(beta1 = 100, beta2 = 1, beta3 = 2),
                        control = nls.control(maxiter = 100, tol = 1e-6))},
           'm14' = {nls(y ~ c * (1-exp(-(x / a) ^ b)),
                        start = list(a = 5, b = 1, c = 100),
                        control = nls.control(maxiter = 100, tol = 1e-6))},
           'm16' = {nls(y ~ c * (x / 204.16) ^ b ,
                        start = list( b = 0.26, c = 100),
                        control = nls.control(maxiter = 1000, tol = 10))}
           
           
    )

#plot(log(x),y)
## Modelo ajustado sin efecto aleatorio
# f<-function(beta1,beta2,beta3,x){beta1 * exp(-beta2 * exp(-beta3 * x))}
summary(modelo_inicial)
predito = predict(modelo_inicial,x)

#plot(predito)
# Visualization of initial model
xyplot(y ~ x | ID,  groups=ID,data=dados, main='M13 (sin efectos aleatorios)',xlab = "Size", ylab = "Acumulado Pasante", pch = 16, 
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...);  
           y.prev <- predito[dados$ID == unique(dados$ID)[panel.number()]]
           llines(x, y.prev, lwd = 2, col = 1);  
           panel.text(x = median(x)+54.5, y = median(y), labels = unique(dados$ID)[panel.number()], pos = 1, cex = 1.2);
       })   

#
#xyplot(y ~ x ,data=dados, main='Bertalanffy (sin efectos aleatorios)',xlab = "Size", ylab = "Acumulado Pasante", pch = 16, 
#panel = function(x, y, ...) {
#    panel.xyplot(x, y, ...);  
#    y.prev <- predito
#    llines(x, y.prev, lwd = 2, col = 1);  
#    panel.text(x = median(x)+54.5, y = median(y), pos = 1, cex = 1.2);
#})   

beta_inicial = coef(modelo_inicial) 


#beta_inicial = c(a=204.16,coef(modelo_inicial))
################## las que no corren ################################
#beta_inicial = c(a=4.1733754, b=0.8344896, c=0.5364471) # m9 skaggs
#beta_inicial = c(1.2516313, 0.1349597, 10.5807787) # m11 Richards
#beta_inicial = c(6.760509, 11.187821) # m12 Swebrec
#beta_inicial = c(103.01717680, 0.05593692, 0.01226751) # m13 bass
#beta_inicial = c(16.123846, 1.081139) # m14 Rossin Rammler
#beta_inicial = c(0.1181964, 8.9358750) # m15 gompertz
#beta_inicial = c(36.9306100, 0.1744041) # m16 Gaudin Schuhmann
#beta_inicial = c(20.713963, 2.550739, 1.17870) # m17 fredlund
#####################################################################

AIC(modelo_inicial);BIC(modelo_inicial)

beta1 <- switch(funcao,
                'logistica'=beta_inicial[1],
                'Gompertz'=beta_inicial[1],
                'm2'=beta_inicial[1],
                'm3'=beta_inicial[1],
                'm4'=beta_inicial[1],
                'm5'=beta_inicial[1],
                'm6'=beta_inicial[1],
                'm7'=beta_inicial[1],
                'm8'=beta_inicial[1],
                'm9'=beta_inicial[1],
                'm10'=beta_inicial[1],
                'm11'=beta_inicial[1],
                'm12'=beta_inicial[1],
                'm13'=beta_inicial[1],
                'm14'=beta_inicial[1],
                'm15'=beta_inicial[1],
                'm16'=beta_inicial[1],
                'm17'=beta_inicial[1])

beta2 <- switch(funcao,
                'logistica'=beta_inicial[3],
                'Gompertz'=beta_inicial[2],
                'm2'=beta_inicial[2],
                'm3'=beta_inicial[2],
                'm4'=beta_inicial[2],
                'm5'=beta_inicial[2],
                'm6'=beta_inicial[2],
                'm7'=beta_inicial[2],
                'm8'=beta_inicial[2],
                'm9'=beta_inicial[2],
                'm10'=beta_inicial[2],
                'm11'=beta_inicial[2],
                'm12'=beta_inicial[2],
                'm13'=beta_inicial[2],
                'm14'=beta_inicial[2],
                'm15'=beta_inicial[2],
                'm16'=beta_inicial[2],
                'm17'=beta_inicial[2])

beta3 <- switch(funcao,
                'logistica'=1/beta_inicial[2],
                'Gompertz'=beta_inicial[3],
                'm2'=beta_inicial[3],
                'm3'=beta_inicial[3],
                'm4'=beta_inicial[3],
                'm5'=beta_inicial[3],
                'm6'=beta_inicial[3],
                'm7'=beta_inicial[3],
                'm8'=beta_inicial[3],
                'm9'=beta_inicial[3],
                'm10'=beta_inicial[3],
                'm11'=beta_inicial[3],
                'm12'=beta_inicial[3],
                'm13'=beta_inicial[3],
                'm14'=beta_inicial[3],
                'm15'=beta_inicial[3],
                'm16'=beta_inicial[3],
                'm17'=beta_inicial[3])

beta = c(beta1,beta2,beta3)

## Ajuste de las funciones aleatorias
z<-cbind(rep(1,n),x)  
u<-rep(1,length(ii))  
n = length(y)
p<-3
sigma2<-10;  
d11<-1; d12<-0;  d22<-1;   
D<-matrix(c(d11,d12,d12,d22),nrow=2,byrow=T)  #valores iniciais
xk<-c(sigma2,d11,d12,d22); 
Sigma<-matrix(0,nrow=n,ncol=n) 
ff = f(beta,x,funcao)
JJ = J(beta,x,funcao)

# Calculo de distnacia de mahalanobis
for(i in ii){ Sigma[s(i),s(i)] <- z[s(i),] %*% D %*% t(z[s(i),])+ sigma2 * Imi(m[i]) }
SigmaInv<-qr.solve(Sigma,tol=1e-100000)

for(i in ii){ u[i] <- t((y-ff)[s(i)]) %*% SigmaInv[s(i),s(i)] %*% (y-ff)[s(i)] }

beta<-c(beta1,beta2,beta3)
vetorPar = c(c(beta),sigma2,d11,d12,d22)

atualiza_modelo <- updatemodel(vetorPar,y,ff,modelo)
Wg = atualiza_modelo$Wg
WgLinha = atualiza_modelo$WgLinha
xi = atualiza_modelo$xy
dg = atualiza_modelo$dg
fg = atualiza_modelo$fg


# Elliptical mixed nonlinear model
contador<-0  

ff = f(beta,x,funcao)
JJ = J(beta,x,funcao)

#  r<-y-f        

# Iterative procedure for estimating betas, sigma2, D 
# Fisher-scoring method
controleWhile = 10
while(abs(controleWhile)>1e-8){
    
    contador<-contador+1
    vetorPar0 = vetorPar          
    
    # Fixed-effects estimation
    incr1<-incr2<-rep(0,length(beta))
    
    for(i in 1:length(m)){
        incr1<- incr1 + 4*dg[i]/m[i]*(t(JJ[s(i),])  %*% SigmaInv[s(i),s(i)] %*% JJ[s(i),])  # alpha
        incr2<- incr2 + (-2 * Wg[i]) * t(JJ[s(i),]) %*%  SigmaInv[s(i),s(i)] %*% c(y-ff)[s(i)] # U alpha
    }
    
    beta <-beta + qr.solve(incr1,tol=1e-1000) %*% incr2
    beta1<-beta[1]; beta2<-beta[2]; beta3<-beta[3]
    
    ff = f(beta,x,funcao)
    JJ = J(beta,x,funcao)
    
    vetorPar = c(c(beta),sigma2,d11,d12,d22)
    
    # Variance components estimation  
    xkmais<- xk + qr.solve(MatrizInformacao(vetorPar, modelo), tol=1e-1000) %*% Score(vetorPar,modelo)      
    xk<-xkmais   
    
    sigma2<-xk[1]; d11<-xk[2]; d12<-xk[3]; d22<-xk[4]  
    
    D<-matrix(c(d11,d12,d12,d22),nrow=ncol(z),byrow=T)  
    
    Sigma<-matrix(0,nrow=n,ncol=n)   
    for(i in 1:(length(m))){ Sigma[s(i),s(i)]<- z[s(i),] %*% D %*% t(z[s(i),])+ sigma2 * Imi(m[i]) }
    SigmaInv<-qr.solve(Sigma,tol=1e-1000000)
    
    vetorPar = c(c(beta),sigma2,d11,d12,d22) 
    controleWhile = max(abs(1-vetorPar[which(vetorPar0!=0)]/vetorPar0[which(vetorPar0!=0)]))
    
    atualiza_modelo <- updatemodel(vetorPar,y,ff,modelo)
    Wg = atualiza_modelo$Wg
    WgLinha = atualiza_modelo$WgLinha
    xi = atualiza_modelo$xy
    dg = atualiza_modelo$dg
    fg = atualiza_modelo$fg
    
    #  print(LogVerosElliptical(vetorPar, modelo))
    
}
print(beta)  
print(round(xk,8))
print(t(t(Score(vetorPar,modelo))))   
print(LogVerosElliptical(vetorPar, modelo))

varBeta = qr.solve(incr1, tol=1e-100000)
erroBeta = sqrt(diag(varBeta))                             
zBeta = beta/erroBeta
pBeta = 2*pnorm(-abs(zBeta))   

varSigma = qr.solve(MatrizInformacao(vetorPar, modelo),tol=1e-10)

erroSigma = sqrt(diag(varSigma))                             
zSigma = c(sigma2,d11,d12,d22)/erroSigma
pSigma = 2*pnorm(-abs(zSigma))  

mAux = matrix(0,2,2)

matrizSaida = matrix(0,7,5)
matrizSaida[1:3,1] = c(t(beta))   
matrizSaida[1:3,2] = c(t(erroBeta))
matrizSaida[1:3, 3] = zBeta    
matrizSaida[1:3, 4] = beta-2*erroBeta     
matrizSaida[1:3, 5] = beta+2*erroBeta

matrizSaida[4:7,1] = c(sigma2,d11,d12,d22)
matrizSaida[4:7,2] = c(t(erroSigma))
matrizSaida[4:7, 3] = zSigma
matrizSaida[4:7, 4] = c(sigma2,d11,d12,d22)-2*erroSigma
matrizSaida[4:7, 5] = c(sigma2,d11,d12,d22)+2*erroSigma

dimnames(matrizSaida) = list(c("beta1", "beta2", "beta3", "sigma2","d11", "d12", "d22"), c("Estimativa","Erro","Z","LI","LS"))
print(round(matrizSaida,6)) 

AIC = -2 * LogVerosElliptical(vetorPar, modelo) + 2 * (length(c(beta)) + 1 + nrow(D)*(nrow(D)+1)/2)
BIC = -2 * LogVerosElliptical(vetorPar, modelo) + 2 * (length(c(beta)) + 1 + nrow(D)*(nrow(D)+1)/2)*log(length(m))

print(AIC);print(BIC)

u<-rep(1,length(m))

for(i in 1:(length(m))) {u[i]<- t(y[s(i)] - ff[s(i)])%*% SigmaInv[s(i),s(i)] %*%(y[s(i)] - ff[s(i)])}


# Random effects prediction
b <- matrix(0,ncol=(length(m)),nrow=2)

for(i in 1:(length(m))){ b[,i] <- D %*% t(z[s(i),]) %*% SigmaInv[s(i),s(i)] %*% (y-ff)[s(i)] }

# Predicting observations
y.prev<-y

for(i in 1:(length(m))) {y.prev[s(i)] <- sigma2*SigmaInv[s(i),s(i)]  %*% (ff)[s(i)] + (Imi(m[i])-sigma2*SigmaInv[s(i),s(i)]) %*% (y)[s(i)] }

subject<-kronecker(1:50,rep(1,15))



# Fitted models visualization      
setTrellis()   

xyplot(y.prev ~ x| ID, groups= ID, data=dados, background="white" ,
       type="o", pch=16, cex=0.5, col=1,  xlab="Size", ylab="Acumulado Pasante",  label.curves=FALSE )

# xx <- seq(0, max(x), length.out = 100)
# ff <- beta1/(1+exp(-(xx-beta2)/beta3))

xyplot(y ~ x | ID,  groups=ID, data=dados, main='Modelo Bertalanffy con efectos aleatorios',xlab = "Size", ylab = "Acumulado Pasante", pch = 16, 
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...);  # Plotando os pontos de observação
           y.prev <- y.prev[dados$ID == unique(dados$ID)[panel.number()]]
           #
           #           y.prev2 <- dados$mu[dados$ID == unique(dados$ID)[panel.number()]]
           #
           llines(x, y.prev, lwd = 2, col = 1);  # Plotando as curvas preditas
           panel.text(x = median(x)+50.2, y = median(y), labels = unique(dados$ID)[panel.number()], pos = 1, cex = 1.2);
       })   


###############################
xyplot(y ~ x | ID, groups = ID, data = dados, 
       main = 'Modelo Bertalanffy con efectos aleatorios',
       xlab = "Size", ylab = "Acumulado Pasante", pch = 16, 
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...);  # Plotando los puntos de observación
           y.prev <- y.prev[dados$ID == unique(dados$ID)[panel.number()]]
           #           y.prev2 <- dados$mu[dados$ID == unique(dados$ID)[panel.number()]]
           
           # Plotando las curvas preditas
           llines(x, y.prev, lwd = 2, col = 1);  
           # Añadiendo la línea gris para y.prev2
           #           llines(x, y.prev2, lwd = 2, col = "gray");  # Línea gris para y.prev2
           
           panel.text(x = median(x) + 50.2, y = median(y), 
                      labels = unique(dados$ID)[panel.number()], 
                      pos = 1, cex = 1.2);
       })   
###############################

Sigma<-matrix(0,nrow=n,ncol=n)   
for (i in 1:(length(m))) Sigma[s(i),s(i)]<- z[s(i),] %*% D %*% t(z[s(i),])+ sigma2 * Imi(m[i])                   
SigmaInv<-qr.solve(Sigma,tol=1e-1000000) 
ff = f(beta,x,funcao)
JJ = J(beta,x,funcao)

# Mahalanobis distance

u<-rep(1,(length(m)))
for (i in 1:(length(m))) u[i]<-
    switch(modelo,
           'normal'= t(y[s(i)] - ff[s(i)])%*% SigmaInv[s(i),s(i)] %*%(y[s(i)] - ff[s(i)]),
           't-Student' = t(y[s(i)] - ff[s(i)])%*% SigmaInv[s(i),s(i)] %*%(y[s(i)] - ff[s(i)]) / m[i]
    )

plot(u, pch=16, ylab="Dist. Mahalanobis", xlab= "ID Sample", )    
points(5, 119, col = "blue", pch = 1, bg = "blue", cex = 3, lwd=3)  
points(1, 138, col = "blue", pch = 1, bg = "blue", cex = 3, lwd =3)
points(5, 119, col = "black", pch = 2, bg = "blue", cex = 1, lwd=3)  
points(1, 138, col = "black", pch = 22, bg = "blue", cex = 1, lwd =3)
axis(1, at = 1:15)  # Eje X del 1 al 15


#  identify(u)




##################


residuales <- c((y[s(1)] - ff[s(1)])%*%SigmaInv[s(1),s(1)],
                (y[s(2)] - ff[s(2)])%*%SigmaInv[s(2),s(2)],
                (y[s(3)] - ff[s(3)])%*%SigmaInv[s(3),s(3)],
                (y[s(4)] - ff[s(4)])%*%SigmaInv[s(4),s(4)],
                (y[s(5)] - ff[s(5)])%*%SigmaInv[s(5),s(5)],
                (y[s(6)] - ff[s(6)])%*%SigmaInv[s(6),s(6)],
                (y[s(7)] - ff[s(7)])%*%SigmaInv[s(7),s(7)],
                (y[s(8)] - ff[s(8)])%*%SigmaInv[s(8),s(8)],
                (y[s(9)] - ff[s(9)])%*%SigmaInv[s(9),s(9)],
                (y[s(10)] - ff[s(10)])%*%SigmaInv[s(10),s(10)],
                (y[s(11)] - ff[s(11)])%*%SigmaInv[s(11),s(11)],
                (y[s(12)] - ff[s(12)])%*%SigmaInv[s(12),s(12)],
                (y[s(13)] - ff[s(13)])%*%SigmaInv[s(13),s(13)],
                (y[s(14)] - ff[s(14)])%*%SigmaInv[s(14),s(14)],
                (y[s(15)] - ff[s(15)])%*%SigmaInv[s(15),s(15)])


hist(residuales)

plot(residuales, ylim = c(-2,2))



(y[s(i)] - ff[s(i)])%*%SigmaInv[s(i),s(i)]


# LONG SET to WIDE
ss <- dados
ss$residuales <-residuales

library(data.table)
setDT(ss)
DB1=dcast(ss, X ~ ID, value.var = "residuales")[,-1]

##--PARAMETRIC SPC--
require(qcc)
qcc.options(bg.margin = "white")

# Shewhart's (X-Bar and Range)
XbarChart = qcc(DB1, type = "xbar", center=0)
RChart = qcc(DB1, type = "R", center = 1.2, ylim =c(0, 2.5))

#limitsRchart <- print(RChart$limits)

# 0 y 2.25

# CuSUM
CuSumChart <- cusum(DB1)
summary(CuSumChart)
# EWMA
EWMAChart <- ewma(DB1, lambda=0.2, nsigmas=3, ylim=c(-0.20,0.20))
summary(EWMAChart)
#limitsewma <- print(EWMAChart$limits)


# -0.20 a +0.20

head(ss)
library(ggplot2)

ss$ID=as.factor(ss$ID)
p1=ggplot(ss,aes(x=X,y=Y,group=ID))+
    geom_line(aes(color=ID)) + xlab(
        expression(paste("Granularidad (",mu,"m)"))) + 
    ylab("Acumulado Pasante")+ theme(legend.position="none")

p2=ggplot(ss,aes(x=log(X),y=Y,group=ID))+
    geom_line(aes(color=ID)) + xlab(
        expression(paste("Granularidad [log(",mu,"m)]"))) + 
    ylab("")

cowplot::plot_grid(p1,p2)

GGally::ggpairs(Hidroc[,c("SIZE","Acum..Pasante")])
#########################

# VISUALIZATION
require(ggplot2)
p11=ggplot(ss, aes(x=residuales)) + geom_histogram() + coord_flip() + 
    scale_y_reverse()  + theme(
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p12=ggplot(ss,aes(x=as.factor(X),y=residuales)) + 
    geom_violin(alpha=0.6) + xlab("Sample (unit)") + 
    ylab("Diff(Observed Sample, Estimated Model)") +
    stat_summary(fun = "median",col="red")+ 
    stat_summary(fun = "mean", col="blue")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cowplot::plot_grid(p11,p12,rel_heights = c(2,1),
                   rel_widths = c(1, 3))



########### Seccion de inferencia Bayesiana ######################################



# Cargando paquetes
library(nlme)
library(Hmisc)
library(lattice)
library(rstan)
library(loo)
require(shinystan)
require(bayesplot)
library(readxl) # Leer excel
dados <- read.csv('Hidrociclon.csv')
#





plot(dados$X, dados$Y)
dados$x <- log(dados$X)

################## logistico ############################

######################## Sin efecto aleatorio

scode2 = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 / (1 + exp(beta0 * (beta1 - x)));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

fitlogistica <- stan(model_code = scode2, 
                     data = list(N=703, x=dados[,6], y=dados[,5]),
                     chain=4, iter=2000, verbose = FALSE)

print(fitlogistica)

elogistica <- extract(fitlogistica, permuted = FALSE)


log_lik_logistica=extract_log_lik(fitlogistica, merge_chains = FALSE)

loo_logistica <- loo(log_lik_logistica)
print(loo_logistica)

waic(log_lik_logistica)

######################### Con efecto aleatorio
scode = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 / (1 + exp( (beta0 + b0[id[i]]) * ((beta1+b1[id[i]])-x[i])));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

head(dados)

fitlogistica2 <- stan(model_code = scode, 
                      data = list( N=703, x=dados[,6], y=dados[,5], id=as.integer(dados[,1]), 
                                   K=length(unique(dados[,1])) ),
                      chain=4, iter=2000, verbose = FALSE)


print(fitlogistica2)

log_lik_logistica2=extract_log_lik(fitlogistica2, merge_chains = FALSE)
waic(log_lik_logistica2)
loo_logistica2 <- loo(log_lik_logistica2)

print(loo_logistica2)

loo_compare(loo_logistica,loo_logistica2)


launch_shinystan(fit1)

## extract samples as a list of arrays
e2 <- extract(fit1, permuted = FALSE)

traceplot(fit1)


mcmc_trace(fit1)



#################### Gaudin -schuhmann


######################## Sin efecto aleatorio

scodeschuhmann = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 * (x / beta0) ^ beta1;
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitschuhmann <- stan(model_code = scodeschuhmann, 
                     data = list(N=703, x=dados[,4], y=dados[,5]),
                     chain=4, iter=2000, verbose = FALSE)
print(fitschuhmann)


log_lik_schuhmann=extract_log_lik(fitschuhmann, merge_chains = FALSE)

loo_schuhmann <- loo(log_lik_schuhmann)
print(loo_schuhmann)



# Con efecto aleatorio

######################### Con efecto aleatorio
scodeschuhmann2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 * (x[i]/ (beta0 + b0[id[i]]) ) ^ (beta1+b1[id[i]]);
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"


head(dados)


fitschuhmann2 <- stan(model_code = scodeschuhmann2, 
                      data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                   K=length(unique(dados[,1])) ),
                      chain=4, iter=2000, verbose = FALSE)


print(fitschuhmann2)

log_lik_schuhmann2=extract_log_lik(fitschuhmann2, merge_chains = FALSE)

loo_schuhmann2 <- loo(log_lik_schuhmann2)
print(loo_schuhmann2)



############## Rosin-Rammler

######################## Sin efecto aleatorio

scoderosin = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 * (1-exp(-(x / beta0) ^ beta1));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitrosin <- stan(model_code = scoderosin, 
                 data = list(N=703, x=dados[,4], y=dados[,5]),
                 chain=4, iter=2000, verbose = FALSE)
print(fitrosin)

log_lik_rosin=extract_log_lik(fitrosin, merge_chains = FALSE)

loo_rosin <- loo(log_lik_rosin)
print(loo_rosin)


# Con efecto aleatorio

######################### Con efecto aleatorio
scoderosin2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 * (1-exp(-(x[i] / (beta0 + b0[id[i]])) ^ (beta1+b1[id[i]])));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"


head(dados)


fitrosin2 <- stan(model_code = scoderosin2, 
                  data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                               K=length(unique(dados[,1])) ),
                  chain=4, iter=2000, verbose = FALSE)

log_lik_rosin2=extract_log_lik(fitrosin2, merge_chains = FALSE)

loo_rosin2 <- loo(log_lik_rosin2)
print(loo_rosin2)


########################## Weibull

######################## Sin efecto aleatorio

scodeweibull = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 - ((100 - beta0) * exp((-(beta1 * x)^beta2)));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitweibull <- stan(model_code = scodeweibull, 
                   data = list(N=703, x=dados[,4], y=dados[,5]),
                   chain=4, iter=2000, verbose = FALSE)
print(fitweibull)

log_lik_weibull=extract_log_lik(fitweibull, merge_chains = FALSE)

loo_weibull <- loo(log_lik_weibull)
print(loo_weibull)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodeweibull2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 - ((100 - (beta0+ b0[id[i]])) * exp((-((beta1+ b1[id[i]]) * x[i])^(beta2+ b2[id[i]]))));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

head(dados)

fitweibull2 <- stan(model_code = scodeweibull2, 
                    data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                 K=length(unique(dados[,1])) ),
                    chain=4, iter=2000, verbose = FALSE)

log_lik_weibull2=extract_log_lik(fitweibull2, merge_chains = FALSE)

loo_weibull2 <- loo(log_lik_weibull2)
print(loo_weibull2)



########################## Gompertz


######################## Sin efecto aleatorio

scodegompertz = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100*exp(-exp(-(beta0*x)^beta1));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitgompertz <- stan(model_code = scodegompertz, 
                    data = list(N=703, x=dados[,4], y=dados[,5]),
                    chain=4, iter=2000, verbose = FALSE)
print(fitgompertz)

log_lik_gompertz=extract_log_lik(fitgompertz, merge_chains = FALSE)

loo_gompertz <- loo(log_lik_gompertz)
print(loo_gompertz)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodegompertz2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 * exp(-exp(-((beta0 + b0[id[i]]) * x[i]) ^ (beta1 + b1[id[i]])));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

head(dados)

fitgompertz2 <- stan(model_code = scodegompertz2, 
                     data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                  K=length(unique(dados[,1])) ),
                     chain=4, iter=2000, verbose = FALSE)

print(fitgompertz2)


log_lik_gompertz2=extract_log_lik(fitgompertz2, merge_chains = FALSE)

loo_gompertz2 <- loo(log_lik_gompertz2)
print(loo_gompertz2)






########################## Hill

######################## Sin efecto aleatorio

scodehill = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = (beta2 * beta0^ beta1)/(beta0^beta1 + x^beta1);
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fithill <- stan(model_code = scodehill, 
                data = list(N=703, x=dados[,4], y=dados[,5]),
                chain=4, iter=2000, verbose = FALSE)
print(fithill)

log_lik_hill=extract_log_lik(fithill, merge_chains = FALSE)

loo_hill <- loo(log_lik_hill)
print(loo_hill)

######################### Con efecto aleatorio
scodehill2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = ((beta2 + b2[id[i]]) * (beta0 + b0[id[i]]) ^ (beta1 + b1[id[i]]))/((beta0 + b0[id[i]]) ^ (beta1 + b1[id[i]]) + x[i] ^ (beta1 + b1[id[i]]));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

head(dados)



fithill2 <- stan(model_code = scodehill2, 
                 data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                              K=length(unique(dados[,1])) ),
                 chain=4, iter=2000, verbose = FALSE)

print(fithill2)


log_lik_hill2=extract_log_lik(fithill2, merge_chains = FALSE)

loo_hill2 <- loo(log_lik_hill2)
print(loo_hill2)


# quede hasta aca

########################## Bass

######################## Sin efecto aleatorio

scodebass = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = beta0 * (1 - exp(-(beta1 + beta2) * x)) ./ (1 + (beta2 / beta1) * exp(-(beta1 + beta2) * x));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitbass <- stan(model_code = scodebass, 
                data = list(N=703, x=dados[,4], y=dados[,5]),
                chain=4, iter=2000, verbose = FALSE)
print(fitbass)

log_lik_bass=extract_log_lik(fitbass, merge_chains = FALSE)

loo_bass <- loo(log_lik_bass)
print(loo_bass)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodebass2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = (beta0 + b0[id[i]]) * (1 - exp(-((beta1 + b1[id[i]]) + (beta2 + b2[id[i]])) * x[i])) ./ (1 + ((beta2 + b2[id[i]]) / (beta1 + b1[id[i]])) * exp(-((beta1 + b1[id[i]]) + (beta2 + b2[id[i]])) * x[i]));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"



head(dados)

fitbass2 <- stan(model_code = scodebass2, 
                 data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                              K=length(unique(dados[,1])) ),
                 chain=4, iter=2000, verbose = FALSE)

print(fitbass2)


log_lik_bass2=extract_log_lik(fitbass2, merge_chains = FALSE)

loo_bass2 <- loo(log_lik_bass2)
print(loo_bass2)




########################## Brody

######################## Sin efecto aleatorio

scodebrody = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = beta0 * (1 - ((beta0 - beta1) / beta0) * exp(-beta2 * x));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitbrody <- stan(model_code = scodebrody, 
                 data = list(N=703, x=dados[,4], y=dados[,5]),
                 chain=4, iter=2000, verbose = FALSE)
print(fitbrody)

log_lik_brody=extract_log_lik(fitbrody, merge_chains = FALSE)

loo_brody <- loo(log_lik_brody)
print(loo_brody)


# Con efecto aleatorio

######################### Con efecto aleatorio
scodebrody2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = (beta0 + b0[id[i]]) * (1 - (((beta0 + b0[id[i]]) - (beta1 + b1[id[i]])) / (beta0 + b0[id[i]])) * exp(-(beta2 + b2[id[i]]) * x[i]));;
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"


head(dados)

fitbrody2 <- stan(model_code = scodebrody2, 
                  data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                               K=length(unique(dados[,1])) ),
                  chain=4, iter=2000, verbose = FALSE)

print(fitbrody2)


log_lik_brody2=extract_log_lik(fitbrody2, merge_chains = FALSE)

loo_brody2 <- loo(log_lik_brody2)
print(loo_brody2)




########################## Bertalanffy

######################## Sin efecto aleatorio

scodebertalanffy = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 * (1 - exp(-beta0 * (x - beta1)));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitbertalanffy <- stan(model_code = scodebertalanffy, 
                       data = list(N=703, x=dados[,4], y=dados[,5]),
                       chain=4, iter=2000, verbose = FALSE)
print(fitbertalanffy)

log_lik_bertalanffy=extract_log_lik(fitbertalanffy, merge_chains = FALSE)

loo_bertalanffy <- loo(log_lik_bertalanffy)
print(loo_bertalanffy)


# Con efecto aleatorio

######################### Con efecto aleatorio
scodebertalanffy2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 * (1 - exp(-(beta0 + b0[id[i]]) * (x[i] - (beta1 + b1[id[i]]))));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

fitbertalanffy2 <- stan(model_code = scodebertalanffy2, 
                        data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                     K=length(unique(dados[,1])) ),
                        chain=4, iter=2000, verbose = FALSE)

print(fitbertalanffy2)


log_lik_bertalanffy2=extract_log_lik(fitbertalanffy2, merge_chains = FALSE)

loo_bertalanffy2 <- loo(log_lik_bertalanffy2)
print(loo_bertalanffy2)



########################## Monomolecular

######################## Sin efecto aleatorio

scodemonomolecular = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = beta0 * (1 - ((beta0 - beta1) / beta0) * exp(-beta2 * x));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitmonomolecular <- stan(model_code = scodemonomolecular, 
                         data = list(N=703, x=dados[,4], y=dados[,5]),
                         chain=4, iter=2000, verbose = FALSE)
print(fitmonomolecular)

log_lik_monomolecular=extract_log_lik(fitmonomolecular, merge_chains = FALSE)

loo_monomolecular <- loo(log_lik_monomolecular)
print(loo_monomolecular)


# Con efecto aleatorio

######################### Con efecto aleatorio
scodemonomolecular2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = (beta0 + b0[id[i]]) * (1 - (((beta0 + b0[id[i]]) - (beta1 + b1[id[i]])) / (beta0 + b0[id[i]])) * exp(-(beta2 + b2[id[i]]) * x[i]));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"



fitmonomolecular2 <- stan(model_code = scodemonomolecular2, 
                          data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                       K=length(unique(dados[,1])) ),
                          chain=4, iter=2000, verbose = FALSE)

print(fitmonomolecular2)


log_lik_monomolecular2=extract_log_lik(fitmonomolecular2, merge_chains = FALSE)

loo_monomolecular2 <- loo(log_lik_monomolecular2)
print(loo_monomolecular2)


########################## Lima

######################## Sin efecto aleatorio

scodelima = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = beta0 + (100 - beta0) / (1 + (beta1 / x) ^ beta2) ^ (1 - 1 / beta2);  
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitlima <- stan(model_code = scodelima, 
                data = list(N=703, x=dados[,4], y=dados[,5]),
                chain=4, iter=2000, verbose = FALSE)
print(fitlima)

log_lik_lima=extract_log_lik(fitlima, merge_chains = FALSE)

loo_lima <- loo(log_lik_lima)
print(loo_lima)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodelima2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = (beta0 + b0[id[i]]) + (100 - (beta0 + b0[id[i]])) / (1 + ((beta1 + b1[id[i]]) / x[i]) ^ (beta2 + b2[id[i]])) ^ (1 - 1 / (beta2 + b2[id[i]]));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

fitlima2 <- stan(model_code = scodelima2, 
                 data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                              K=length(unique(dados[,1])) ),
                 chain=4, iter=2000, verbose = FALSE)

print(fitlima2)


log_lik_lima2=extract_log_lik(fitlima2, merge_chains = FALSE)

loo_lima2 <- loo(log_lik_lima2)
print(loo_lima2)





########################## Skaggs

######################## Sin efecto aleatorio

scodeskaggs = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 / (1 + (100 / beta0 - 1) * exp(-beta1 * x ^ beta2));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitskaggs <- stan(model_code = scodeskaggs, 
                  data = list(N=703, x=dados[,4], y=dados[,5]),
                  chain=4, iter=2000, verbose = FALSE)
print(fitskaggs)

log_lik_skaggs=extract_log_lik(fitskaggs, merge_chains = FALSE)

loo_skaggs <- loo(log_lik_skaggs)
print(loo_skaggs)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodeskaggs2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 / (1 + (100 / (beta0 + b0[id[i]]) - 1) * exp(-(beta1 + b1[id[i]]) * x[i] ^ (beta2 + b2[id[i]])));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

fitskaggs2 <- stan(model_code = scodeskaggs2, 
                   data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                K=length(unique(dados[,1])) ),
                   chain=4, iter=2000, verbose = FALSE)

print(fitskaggs2)


log_lik_skaggs2=extract_log_lik(fitskaggs2, merge_chains = FALSE)

loo_skaggs2 <- loo(log_lik_skaggs2)
print(loo_skaggs2)



########################## Fredlund

######################## Sin efecto aleatorio

scodefredlund = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  //mu = 100 / (log(exp(1) + ((beta0 / x) ^ beta1) ^ beta2) * (1 - (log(1 + 0.105 / x) / log(1 + 0.105 / 0.00001)) ^ 7));
  mu = 100 ./ (log(exp(1) + pow(beta0 ./ x, beta1 * beta2)) .* (1 - pow(log(1 + 0.105 ./ x) ./ log(1 + 0.105 / 0.00001), 7))); // version rstan
  
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"
min(x)
head(dados)
fitfredlund <- stan(model_code = scodefredlund, 
                    data = list(N=703, x=dados[,4], y=dados[,5]),
                    chain=4, iter=2000, verbose = FALSE)
print(fitfredlund)

log_lik_fredlund=extract_log_lik(fitfredlund, merge_chains = FALSE)

loo_fredlund <- loo(log_lik_fredlund)
print(loo_fredlund)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodefredlund2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 ./ (log(exp(1) + pow((beta0 + b0[id[i]]) ./ x[i], (beta1 + b1[id[i]]) * (beta2 + b2[id[i]]))) .* (1 - pow(log(1 + 0.105 ./ x[i]) ./ log(1 + 0.105 / 0.00001), 7)));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"


fitfredlund2 <- stan(model_code = scodefredlund2, 
                     data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                  K=length(unique(dados[,1])) ),
                     chain=4, iter=2000, verbose = FALSE)

print(fitfredlund2)


log_lik_fredlund2=extract_log_lik(fitfredlund2, merge_chains = FALSE)

loo_fredlund2 <- loo(log_lik_fredlund2)
print(loo_fredlund2)




########################## Havercamp

######################## Sin efecto aleatorio

scodehavercamp = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 / (1 + (beta0 / x) ^ beta1) ^ beta2;
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fithavercamp <- stan(model_code = scodehavercamp, 
                     data = list(N=703, x=dados[,4], y=dados[,5]),
                     chain=4, iter=2000, verbose = FALSE)
print(fithavercamp)

log_lik_havercamp=extract_log_lik(fithavercamp, merge_chains = FALSE)

loo_havercamp <- loo(log_lik_havercamp)
print(loo_havercamp)



# Con efecto aleatorio

######################### Con efecto aleatorio
scodehavercamp2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 / (1 + ((beta0 + b0[id[i]]) / x[i]) ^ (beta1 + b1[id[i]])) ^ (beta2 + b2[id[i]]);
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

fithavercamp2 <- stan(model_code = scodehavercamp2, 
                      data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                   K=length(unique(dados[,1])) ),
                      chain=4, iter=2000, verbose = FALSE)



print(fithavercamp2)


log_lik_havercamp2=extract_log_lik(fithavercamp2, merge_chains = FALSE)

loo_havercamp2 <- loo(log_lik_havercamp2)
print(loo_havercamp2)





########################## Morganetal

######################## Sin efecto aleatorio

scodemorganetal = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 - ((100 - beta0) / (1 + (beta1 * x) ^ beta2));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitmorganetal <- stan(model_code = scodemorganetal, 
                      data = list(N=703, x=dados[,4], y=dados[,5]),
                      chain=4, iter=2000, verbose = FALSE)
print(fitmorganetal)

log_lik_morganetal=extract_log_lik(fitmorganetal, merge_chains = FALSE)

loo_morganetal <- loo(log_lik_morganetal)
print(loo_morganetal)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodemorganetal2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 - ((100 - (beta0 + b0[id[i]])) / (1 + ((beta1 + b1[id[i]]) * x[i]) ^ (beta2 + b2[id[i]])));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

fitmorganetal2 <- stan(model_code = scodemorganetal2, 
                       data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                    K=length(unique(dados[,1])) ),
                       chain=4, iter=2000, verbose = FALSE)

print(fitmorganetal2)


log_lik_morganetal2=extract_log_lik(fitmorganetal2, merge_chains = FALSE)

loo_morganetal2 <- loo(log_lik_morganetal2)
print(loo_morganetal2)




########################## Richards

######################## Sin efecto aleatorio

scoderichards = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> beta2;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 * (1+(beta0 - 1) * exp(- beta1 *(x-beta2)))^(1/(1-beta0));
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
  beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitrichards <- stan(model_code = scoderichards, 
                    data = list(N=703, x=dados[,4], y=dados[,5]),
                    chain=4, iter=2000, verbose = FALSE)
print(fitrichards)

log_lik_richards=extract_log_lik(fitrichards, merge_chains = FALSE)

loo_richards <- loo(log_lik_richards)
print(loo_richards)


# Con efecto aleatorio

######################### Con efecto aleatorio
scoderichards2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    real beta2;
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    vector[K] b2;
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 * (1+((beta0+ b0[id[i]]) - 1) * exp(- (beta1+ b1[id[i]]) *(x[i]-(beta2+ b2[id[i]]))))^(1/(1-(beta0+ b0[id[i]])));
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
    beta2 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"

fitrichards2 <- stan(model_code = scoderichards2, 
                     data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                  K=length(unique(dados[,1])) ),
                     chain=4, iter=2000, verbose = FALSE)

print(fitrichards2)


log_lik_richards2=extract_log_lik(fitrichards2, merge_chains = FALSE)

loo_richards2 <- loo(log_lik_richards2)
print(loo_richards2)






########################## Swebrec
######################## Sin efecto aleatorio

scodeswebrec = "// BLOQUE DE LOS DATOS
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
// BLOQUE DE LOS PARAMETROS
parameters {
  real<lower=0> beta0;
  real<lower=0> beta1;
  real<lower=0> sigma;
}
// BLOQUE DE LA TRANSFORMADA
transformed parameters {
  vector[N] mu;
  mu = 100 / (1+(log(max(x)/x)/log(max(x))/beta1))^beta0;
    
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
  //VEROSIMILITUD
  y ~ normal(mu, sigma);
  //PRIORI
  beta0 ~ normal(1, 100);
  beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}
"

head(dados)
fitswebrec <- stan(model_code = scodeswebrec, 
                   data = list(N=703, x=dados[,4], y=dados[,5]),
                   chain=4, iter=2000, verbose = FALSE)
print(fitswebrec)

log_lik_swebrec=extract_log_lik(fitswebrec, merge_chains = FALSE)

loo_swebrec <- loo(log_lik_swebrec)
print(loo_swebrec)

# Con efecto aleatorio

######################### Con efecto aleatorio
scodeswebrec2 = "// BLOQUE DE LOS DATOS
data {
    int<lower=1> N;
    vector[N] x;
    vector[N] y;
    int<lower=1> K;
    int<lower=1, upper=K> id[N];
}
// BLOQUE DE LOS PARAMETROS
parameters {
    real beta0;
    real beta1;      // or vector[2] beta; //intercept and slope
    vector[K] b0;
    vector[K] b1;    // or matrix[2,K] w; //groups intercepts, slopes
    real<lower=0> sigma;
}
transformed parameters {
    vector[N] mu;
    for (i in 1:N){
    mu[i] = 100 / (1+(log(max(x)/x[i])/log(max(x))/(beta1+b1[id[i]])))^(beta0 + b0[id[i]]);
    }
}
// BLOQUE MODELO (PRIORI + VEROSIMILITUD)
model {
    //VEROSIMILITUD
    y ~ normal(mu, sigma);
    //PRIORI
    beta0 ~ normal(1, 100);
    beta1 ~ normal(0.2, 100);
}
// GENERATED QUATITIES BLOCK
generated quantities {
    array[N] real log_lik;
    for (n in 1:N) {
        log_lik[n] = normal_rng(mu[n],sigma);
    }
}

"


fitswebrec2 <- stan(model_code = scodeswebrec2, 
                    data = list( N=703, x=dados[,4], y=dados[,5], id=as.integer(dados[,1]), 
                                 K=length(unique(dados[,1])) ),
                    chain=4, iter=2000, verbose = FALSE)

print(fitswebrec2)


log_lik_swebrec2=extract_log_lik(fitswebrec2, merge_chains = FALSE)

loo_swebrec2 <- loo(log_lik_swebrec2)
print(loo_swebrec2)







####################################################


loo_compare(loo_logistica,loo_schuhmann,loo_rosin,loo_weibull,
            loo_gompertz,loo_hill,loo_bass,loo_brody,loo_bertalanffy
            ,loo_monomolecular,loo_lima,loo_skaggs,loo_fredlund
            ,loo_havercamp,loo_morganetal,loo_richards,loo_swebrec)

# bertalanffy
# lima
# rosin
# weibull
# skaggs
# havecamp


loo_compare(loo_logistica2,loo_schuhmann2,loo_rosin2,loo_weibull2,
            loo_gompertz2,loo_hill2,loo_bass2,loo_brody2,loo_bertalanffy2
            ,loo_monomolecular2,loo_lima2,loo_skaggs2,loo_fredlund2
            ,loo_havercamp2,loo_morganetal2,loo_richards2,loo_swebrec2)

# TODOS (TOP 5)
loo_compare(loo_logistica, # 1
            loo_schuhmann, # 2
            loo_rosin, # 3
            loo_weibull, # 4
            loo_gompertz, # 5
            loo_hill, # 6
            loo_bass, # 7
            loo_brody, #8
            loo_bertalanffy, # 9
            loo_monomolecular, # 10
            loo_lima, # 11
            loo_skaggs, # 12
            loo_fredlund, # 13
            loo_havercamp, # 14
            loo_morganetal, # 15
            loo_richards, # 16
            loo_swebrec, # 17
            loo_logistica2, # 18
            loo_schuhmann2, # 19
            loo_rosin2, # 20
            loo_weibull2, # 21
            loo_gompertz2, # 22
            loo_hill2, # 23
            loo_bass2, # 24
            loo_brody2, # 25
            loo_bertalanffy2, # 26
            loo_monomolecular2, # 27
            loo_lima2, # 28
            loo_skaggs2, # 29
            loo_fredlund2, # 30
            loo_havercamp2, # 31
            loo_morganetal2, # 32
            loo_richards2, # 33
            loo_swebrec2) # 34


# Mejor Havercamp2
# hill2
# rosin2
# morgan et al2
# logistica2
# fredlund
# bertalanffy
# lima 
# skaggs
# havercamp


# el modelo 17 es el logistico

loo_compare(
    waic(log_lik_logistica),# 1
    waic(log_lik_logistica2),# 2
    waic(log_lik_bass),# 3
    waic(log_lik_bass2),# 4
    waic(log_lik_bertalanffy), # 5
    waic(log_lik_bertalanffy2), # 6
    waic(log_lik_brody), # 7
    waic(log_lik_brody2), # 8
    waic(log_lik_fredlund), # 9
    waic(log_lik_fredlund2), # 10
    waic(log_lik_gompertz), # 11
    waic(log_lik_gompertz2), # 12
    waic(log_lik_havercamp), # 13
    waic(log_lik_havercamp2), # 14
    waic(log_lik_hill), # 15
    waic(log_lik_hill2), # 16
    waic(log_lik_lima), # 17
    waic(log_lik_lima2), # 18
    waic(log_lik_monomolecular), # 19
    waic(log_lik_monomolecular2), # 20
    waic(log_lik_morganetal), # 21
    waic(log_lik_morganetal2), # 22
    waic(log_lik_richards), # 23
    waic(log_lik_richards2), # 24
    waic(log_lik_rosin), # 25
    waic(log_lik_rosin2), # 26
    waic(log_lik_schuhmann), # 27
    waic(log_lik_schuhmann2), # 28
    waic(log_lik_skaggs), # 29
    waic(log_lik_skaggs2), # 30
    waic(log_lik_swebrec), # 31
    waic(log_lik_swebrec2), # 32
    waic(log_lik_weibull), # 33
    waic(log_lik_weibull2) # 34
)


#### WAIC

loo(log_lik_logistica)$looic; waic(log_lik_logistica)$waic
loo(log_lik_logistica2)$looic; waic(log_lik_logistica2)$waic


loo(log_lik_bass)$looic; waic(log_lik_bass)$waic
loo(log_lik_bass2)$looic; waic(log_lik_bass2)$waic

loo(log_lik_bertalanffy)$looic; waic(log_lik_bertalanffy)$waic
loo(log_lik_bertalanffy2)$looic; waic(log_lik_bertalanffy2)$waic


loo(log_lik_brody)$looic; waic(log_lik_brody)$waic
loo(log_lik_brody2)$looic; waic(log_lik_brody2)$waic

loo(log_lik_fredlund)$looic; waic(log_lik_fredlund)$waic
loo(log_lik_fredlund2)$looic; waic(log_lik_fredlund2)$waic

loo(log_lik_gompertz)$looic; waic(log_lik_gompertz)$waic
loo(log_lik_gompertz2)$looic; waic(log_lik_gompertz2)$waic

loo(log_lik_havercamp)$looic; waic(log_lik_havercamp)$waic
loo(log_lik_havercamp2)$looic; waic(log_lik_havercamp2)$waic

loo(log_lik_hill)$looic; waic(log_lik_hill)$waic
loo(log_lik_hill2)$looic; waic(log_lik_hill2)$waic

loo(log_lik_lima)$looic; waic(log_lik_lima)$waic
loo(log_lik_lima2)$looic; waic(log_lik_lima2)$waic

loo(log_lik_monomolecular)$looic; waic(log_lik_monomolecular)$waic
loo(log_lik_monomolecular2)$looic; waic(log_lik_monomolecular2)$waic

loo(log_lik_morganetal)$looic; waic(log_lik_morganetal)$waic
loo(log_lik_morganetal2)$looic; waic(log_lik_morganetal2)$waic

loo(log_lik_richards)$looic; waic(log_lik_richards)$waic
loo(log_lik_richards2)$looic; waic(log_lik_richards2)$waic

loo(log_lik_rosin)$looic; waic(log_lik_rosin)$waic
loo(log_lik_rosin2)$looic; waic(log_lik_rosin2)$waic


loo(log_lik_schuhmann)$looic; waic(log_lik_schuhmann)$waic
loo(log_lik_schuhmann2)$looic; waic(log_lik_schuhmann2)$waic

loo(log_lik_skaggs)$looic; waic(log_lik_skaggs)$waic
loo(log_lik_skaggs2)$looic; waic(log_lik_skaggs2)$waic

loo(log_lik_swebrec)$looic; waic(log_lik_swebrec)$waic
loo(log_lik_swebrec2)$looic; waic(log_lik_swebrec2)$waic

loo(log_lik_weibull)$looic; waic(log_lik_weibull)$waic
loo(log_lik_weibull2)$looic; waic(log_lik_weibull2)$waic



library(rstanarm)
library(plotly)
library(dplyr)
library(tidyverse)
library(lattice)
library(rstan)
library(dplyr)


b <- extract(fithavercamp2)
b <- rstan::extract(fithavercamp2)


vector_mu <- as.vector(b$mu)

promedios <- tapply(vector_mu, (seq_along(vector_mu) - 1) %/% 4000, mean)
promedios <- as.vector(promedios)
length(promedios)
dados <- dados %>%
    mutate(mu = promedios)


