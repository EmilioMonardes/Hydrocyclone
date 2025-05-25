updatemodel<-function(vetorPar,y,f,modelo){
  switch(modelo,
         'normal' = {
           mm <- m
           dg <- mm / 4
           fg <- mm * (mm + 2) / 4
           xi <- rep(1,length(m))
           Wg <- rep(-0.5,length(m))
           WgLinha <- rep(0,length(m))
           return(list(Wg = Wg, WgLinha=WgLinha, xi = xi, dg = dg, fg = fg))
         },
         't-Student' = {
           for (i in 1:length(m)) u[i] <- t((y-ff)[s(i)]) %*% SigmaInv[s(i),s(i)] %*% (y-ff)[s(i)]                                
           Wg = -(nu+m)/(2*(nu+u)) 
           WgLinha = (nu+m)/(2*(nu+u)^2) 
           xi = nu/(nu-2)
           dg = (m/4)*((nu+m)/(nu+m+2))
           fg = (m*(m+2)/4)*((nu+m)/(nu+m+2)) 
           return(list(Wg = Wg, WgLinha=WgLinha, xi = xi, dg = dg, fg = fg))
         },
         stop("Modelo não reconhecido")
  )
}


LogVerosElliptical = function(t, modelo)
{
  beta<-t[1:6]; sigma2<-t[4]; d11<-t[5]; d12<-t[6]; d22<-t[7] 
  beta1<-beta[1]; beta2<-beta[2]; beta3<-100

  ff = f(beta,x,funcao)
  JJ = J(beta,x,funcao)
   
  z<-cbind(rep(1,n),x)

  D<-matrix(c(d11,d12,d12,d22),nrow=ncol(z),byrow=T)
  
  atualiza_modelo <- updatemodel(vetorPar,y,f,modelo)
  Wg = atualiza_modelo$Wg
  WgLinha = atualiza_modelo$WgLinha
  xi = atualiza_modelo$xy
  dg = atualiza_modelo$dg
  fg = atualiza_modelo$fg
  
  Sigma<-matrix(0,nrow=n,ncol=n)     
  c1<-c2<-rep(0,(length(m)))   
  
  for (i in 1:(length(m))) Sigma[s(i),s(i)] <- z[s(i),] %*% D %*% t(z[s(i),])+ sigma2 * Imi(m[i]) 
    SigmaInv<-qr.solve(Sigma,tol=1e-100000)                                                                                          
  
  for (i in 1:(length(m))) {
    c1[i] <- (y-ff)[s(i)] %*% SigmaInv[s(i),s(i)] %*% (y-ff)[s(i)]
    c2[i] <- determinant(Sigma[s(i),s(i)])$modulus[1]  
  }   
  
  switch(modelo,
         'normal' = {-((n)/2)*log(2*pi)-1/2*sum(c2)-1/2*sum(c1)},
         't-Student' = {sum(lgamma((m+nu)/2)-lgamma(nu/2)-m/2*log(nu*pi))-1/2*sum(c2)-sum(log((1+1/nu*c1)^(((m+nu)/2))))},
         stop("Modelo não reconhecido")
  )
  
}      

Score = function(t, modelo)
{
  beta1<-t[1]; beta2<-t[2]; beta3<-t[3]
  sigma2<-t[4]; d11<-t[5]; d12<-t[6]; d22<-t[7] 
  
  atualiza_modelo <- updatemodel(vetorPar,y,f,modelo)
  Wg = atualiza_modelo$Wg
  WgLinha = atualiza_modelo$WgLinha
  xi = atualiza_modelo$xy
  dg = atualiza_modelo$dg
  fg = atualiza_modelo$fg
  
  Sigma<-matrix(0,nrow=n,ncol=n)     
  c11<-c12<-c13<-c14<-rep(0,(length(m)))
  for (i in 1:(length(m)))                  
    Sigma[s(i),s(i)]<- z[s(i),] %*% D %*% t(z[s(i),])+ sigma2 * Imi(m[i])                 
  SigmaInv<-qr.solve(Sigma, tol=1e-100000) 
  for (i in 1:(length(m)))
  {   
    c11[i] <- sum(diag(SigmaInv[s(i),s(i)])) + 2*Wg[i] * t((y-ff)[s(i)]) %*% SigmaInv[s(i),s(i)] %*% SigmaInv[s(i),s(i)]  %*% (y-ff)[s(i)]
    c12[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1]))))+ 2*Wg[i] * t((y-ff)[s(i)]) %*% SigmaInv[s(i),s(i)] %*%  (z[s(i),1] %*% t(z[s(i),1])) %*% SigmaInv[s(i),s(i)]  %*% (y-ff)[s(i)] 
    c13[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+(z[s(i),2] %*% t(z[s(i),1])))))+ 2*Wg[i] * t((y-ff)[s(i)]) %*% SigmaInv[s(i),s(i)] %*%  ((z[s(i),1] %*% t(z[s(i),2]))+(z[s(i),2] %*% t(z[s(i),1]))) %*% SigmaInv[s(i),s(i)]  %*% (y-ff)[s(i)]
    c14[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))+ 2*Wg[i] * t((y-ff)[s(i)]) %*% SigmaInv[s(i),s(i)] %*%  (z[s(i),2] %*% t(z[s(i),2])) %*% SigmaInv[s(i),s(i)]  %*% (y-ff)[s(i)]
  }    
  c(-1/2*sum(c11), -1/2*sum(c12),  -1/2*sum(c13),  -1/2*sum(c14))
}    

MatrizInformacao<-function(t,modelo)
{
  beta1<-t[1]; beta2<-t[2]; beta3<-t[3]
  sigma2<-t[4]; d11<-t[5]; d12<-t[6]; d22<-t[7] 
  
  atualiza_modelo <- updatemodel(vetorPar,y,f,modelo)
  Wg = atualiza_modelo$Wg
  WgLinha = atualiza_modelo$WgLinha
  xi = atualiza_modelo$xy
  dg = atualiza_modelo$dg
  fg = atualiza_modelo$fg
  
  D<-matrix(c(d11,d12,d12,d22),nrow=ncol(z),byrow=T)
  Sigma<-matrix(0,nrow=n,ncol=n)     
  c11<-c12<-c13<-c14<-rep(0,(length(m)))    
  brsi11<-brsi12<-brsi13<-brsi14<-brsi22<-brsi23<-brsi24<-brsi33<-brsi34<-brsi44<-rep(0,(length(m)))     
  for (i in 1:(length(m)))
  {
    brsi11[i] <- sum(diag(SigmaInv[s(i),s(i)]))^2
    brsi12[i] <- sum(diag(SigmaInv[s(i),s(i)]))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1]))))
    brsi13[i] <- sum(diag(SigmaInv[s(i),s(i)]))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),2]) + (z[s(i),2] %*% t(z[s(i),1])))))
    brsi14[i] <- sum(diag(SigmaInv[s(i),s(i)]))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))
    brsi22[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1]))))^2
    brsi23[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1])))) * sum(diag(SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1]))))))
    brsi24[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1])))) * sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))
    brsi33[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1]))))))^2
    brsi34[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1])))))) * sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))
    brsi44[i] <- sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))^2
  }
  K11<-K12<-K13<-K14<-K22<-K23<-K24<-K33<-K34<-K44<-U1<-U2<-U3<-U4<-rep(0,(length(m)))
  for (i in 1:(length(m)))
  {
    K11[i] <- brsi11[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% SigmaInv[s(i),s(i)]))
    K12[i] <- brsi12[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1]))))
    K13[i] <- brsi13[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),2]) + (z[s(i),2] %*% t(z[s(i),1])))))
    K14[i] <- brsi14[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))
    K22[i] <- brsi22[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1])) %*% SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1]))))
    K23[i] <- brsi23[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1])) %*% SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1]))))))
    K24[i] <- brsi24[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),1] %*% t(z[s(i),1])) %*% SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))
    K33[i] <- brsi33[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1])))) %*% SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1]))))))
    K34[i] <- brsi34[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2])) %*% SigmaInv[s(i),s(i)] %*% ((z[s(i),1] %*% t(z[s(i),2]))+((z[s(i),2] %*% t(z[s(i),1]))))))
    K44[i] <- brsi44[i]/4*(4*fg[i]/(mm[i]*(mm[i]+2))-1)+2*fg[i]/(mm[i]*(mm[i]+2))*sum(diag(SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2])) %*% SigmaInv[s(i),s(i)] %*% (z[s(i),2] %*% t(z[s(i),2]))))
  } 
  matrix(c(sum(K11), sum(K12), sum(K13), sum(K14), sum(K12), sum(K22), sum(K23), sum(K24), sum(K13), sum(K23), sum(K33), sum(K34), sum(K14), sum(K24), sum(K34), sum(K44)),nrow=4,ncol=4,byrow=T)
}

s<-function(i)
{
  ma[i]:(ma[i+1]-1)
}

tr<-function(A)
{
  sum(diag(A))
}


Imi <- function(i){ diag(rep(1,i))}



