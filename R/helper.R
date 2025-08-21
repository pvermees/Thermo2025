# pop = the population, i.e. an integer from 1 to 6
# plot = either 'none', 'ages' or 'dates'
# method = either 'fissiontracks' or 'U-Th-He'
# tt = a vector of size n with (central) ages
# sigma = a vector of size n with dispersion factors
# prop = a vector of size n with proportions of populations
# trunc = a logical vector of size n indicating wether the 
#         the lower tail of the component should be truncated
# zeta, rhoD = fission track parameters                                        
# N = geometric mean Ns + Ni and its dispersion
# ThU = geometric Th/U ratio and its dispersion
# nU = geometric mean U cps and its dispersion
# Note: tt, sigma and prop are only used if pop does NOT equal 1...6
getpop <- function(pop=1,tt,sigma,prop,
                   plot='none',method='fissiontracks',
                   zeta=350e-6,rhoD=2.5e6,N=c(100,0.3),
                   ThU=c(1,0.1),nU=c(10000,0.3),
                   from=20,to=500){
    if (pop==1){
        if (missing(tt)) tt <- 100
        sigma <- 0
        prop <- 1
        trunc <- FALSE
    } else if (pop==2){
        if (missing(tt)) tt <- 100
        if (missing(sigma)) sigma <- 0.5
        prop <- 1
        trunc <- FALSE
    } else if (pop==3){
        if (missing(tt) || length(tt)<2) tt <- c(100,200)
        if (missing(sigma) || length(sigma)<2) sigma <- c(0,0)
        if (missing(prop) || length(prop)<2) prop <- c(1,1)/2
        trunc <- c(FALSE,FALSE)
    } else if (pop==4) {
        if (missing(tt) || length(tt)<2) tt <- c(50,100)
        if (missing(sigma) || length(sigma)<2) sigma <- c(0,0.5)
        if (missing(prop) || length(prop)<2) prop <- c(0.5,0.5)
        trunc <- c(FALSE,TRUE)
    } else if (pop==5){
        if (missing(tt) || length(tt)<2) tt <- c(100,300)
        if (missing(sigma) || length(sigma)<2) sigma <- c(0.2,0.3)
        if (missing(prop) || length(prop)<2) prop <- c(0.3,0.7)
        trunc <- c(FALSE,FALSE)
    } else if (pop==6){
        if (missing(tt) || length(tt)<2) tt <- c(100,300)
        if (missing(sigma) || length(sigma)<2) sigma <- c(0.2,0.3)
        if (missing(prop) || length(prop)<2) prop <- c(0.7,0.3)
        trunc <- c(FALSE,FALSE)
    }
    prop <- prop/sum(prop) # normalise
    xx <- exp(seq(from=log(from),to=log(to),length.out=200))
    if (pop==4){
        xx <- sort(c(xx,tt[1],tt[1]+1e-5))
    }
    if (plot %in% c('ages','dates')){
        if (identical(method,'fissiontracks')){
            mut <- t2mu(tt=tt,zeta=zeta,rhoD=rhoD)
            mux <- t2mu(tt=xx,zeta=zeta,rhoD=rhoD)
        } else if (identical(method,'U-Th-He')){
            mut <- log(tt)
            mux <- log(xx)
        } else {
            stop('Unrecognised method.')
        }
        if (pop %in% c(1,3)){
            yy <- seq(from=0,to=1,length.out=length(xx))
        } else if (pop==2){
            yy <- dnorm(mux,mean=mut,sd=sigma)
        } else if (pop==4){
            yy <- dtrunc(log(xx),mean=mut[2],sd=sigma[2],minx=mut[1])
        } else if (pop %in% 5:6){
            np <- length(prop)
            yy <- 0*xx
            for (i in 1:np){
                yy <- yy + prop[i]*dnorm(mux,mean=mut[i],sd=sigma[i])
            }
        }
    }
    if (plot == 'dates'){
        yy <- 0*xx
        for (i in which(sigma>0)){
            pa <- dtrunc(mux,mean=mut[i],sd=sigma[i],
                         minx=ifelse(trunc[i],mut[i-1],-Inf))
            ynew <- 0*yy
            for (j in seq_along(xx)){
                cv <- t2cv(tt=xx[j],method=method,
                           N=N[1],zeta=zeta,rhoD=rhoD,
                           ThU=ThU[1],nU=nU[1])
                ynew <- ynew + pa[j]*dnorm(mux,mean=mux[j],sd=cv)
            }
            yy <- yy + prop[i]*ynew/sum(ynew)
        }
        for (i in which(sigma==0)){
            cv <- t2cv(tt=tt[i],method=method,
                       N=N[1],zeta=zeta,rhoD=rhoD,
                       ThU=ThU[1],nU=nU[1])
            ynew <- dnorm(mux,mean=mut[i],sd=cv)
            yy <- yy + prop[i]*ynew/sum(ynew)
        }
    }
    if (plot=='ages'){
        plot(x=xx,y=yy,type='n',log='x',bty='n',
             xlab='age (Ma)',ylab='',yaxt='n')
        if (pop==1){
            lines(x=c(from,tt,tt,tt,to),y=c(0,0,1,0,0),lwd=2)
        } else if (pop==2){
            lines(x=xx,y=yy,lwd=2)
        } else if (pop==3){
            lines(x=c(from,tt[1],tt[1],tt[1],tt[2],tt[2],tt[2],to),
                  y=c(0,0,1,0,0,1,0,0),lwd=2)
        } else if (pop==4){
            lines(x=xx,y=yy,lwd=2)
            lines(x=c(from,tt[1],tt[1],tt[1],to),
                  y=c(0,0,1,0,0),lwd=2)
        } else if (pop %in% 5:6){
            lines(x=xx,y=yy,lwd=2)
        }
    } else if (plot=='dates'){
        plot(x=xx,y=yy,type='l',log='x',bty='n',
             xlab='date (Ma)',ylab='',yaxt='n',lwd=2)
    }
    invisible(data.frame(t=tt,sigma=sigma,prop=prop,trunc=trunc))
}

# pop = an integer marking a predefined population.
# ng = the desired number of grains
# tt = a vector of size n with (central) ages
# sigma = a vector of size n with dispersion factors
# prop = a vector of size n with proportions of populations
# zeta = two-element vector with the zeta calibration factor
#        and its standard error, in  [yr-1cm2]
# rhoD = two-element vector with the dosimeter track density
#        and its standard error [in cm-2]
# zeta, rhoD = fission track parameters and their standard errors
# N = geometric mean Ns + Ni and its dispersion
# Note: tt, sigma and prop are only used if pop does NOT equal 1...6
FTsamp <- function(pop=0,ng=100,tt=c(50,100),sigma=c(0,0.5),prop=c(0.5,0.5),
                   zeta=c(350e-6,2e-6),rhoD=c(2.5e6,0.05e6),N=c(100,0.3)){
    l38 <- settings('lambda','U238')[1]
    pop <- getpop(pop=pop,tt=tt,sigma=sigma,prop=prop)
    component <- sample(x=1:nrow(pop),size=ng,replace=TRUE)
    mu <- t2mu(tt=pop[component,'t'],zeta=zeta[1],rhoD=rhoD[1])
    truncate <- component %in% which(pop[,'trunc'])
    minx <- rep(-Inf,ng)
    minx[truncate] <- t2mu(tt=pop[component[truncate]-1,'t'],
                           zeta=zeta[1],rhoD=rhoD[1])
    rmu <- rtrunc(ng,mean=mu,sd=pop[component,'sigma'],minx=minx)
    rhosrhoi <- exp(rmu)
    n <- round(exp(rnorm(ng,mean=log(N[1]),sd=N[2])))
    ni <- (n/(1+rhosrhoi))
    ns <- n-ni
    Ni <- rpois(ng,lambda=ni)
    Ns <- rpois(ng,lambda=ns)
    out <- list()
    class(out) <- 'fissiontracks'
    out$x <- cbind(Ns=Ns,Ni=Ni)
    out$rhoD <- rhoD
    out$zeta <- zeta*1e6
    out$format <- 1
    out$t.true <- log(1+0.5*l38*zeta[1]*rhoD[1]*rhosrhoi)/l38
    out
}

# pop = an integer marking a predefined population.
# ng = the desired number of grains
# tt = a vector of size n with (central) ages
# sigma = a vector of size n with dispersion factors
# prop = a vector of size n with proportions of populations
# ThU = geometric mean Th/U ratio and its dispersion
# nU = geometric mean uranium signal (ion counts) and its dispersion
# HeGain = relative sensitivity of noble gass MS relative to ICP-MS
# Note: tt, sigma and prop are only used if pop does NOT equal 1...6
UThHesamp <- function(pop=1,ng=20,tt=c(50,100),sigma=c(0,0.5),prop=c(0.5,0.5),
                      ThU=c(1,0.3),nU=c(500,0.3),HeGain=20){
    pop <- getpop(pop=pop,tt=tt,sigma=sigma,prop=prop)
    component <- sample(x=1:nrow(pop),size=ng,replace=TRUE)
    truncate <- component %in% which(pop[,'trunc'])
    minx <- rep(-Inf,ng)
    minx[truncate] <- log(pop[component[truncate]-1,'t'])
    lt <- rtrunc(ng,mean=log(pop[component,'t']),
                 sd=pop[component,'sigma'],minx=minx)
    tt <- exp(lt)
    lThU <- rnorm(ng,mean=log(ThU[1]),sd=ThU[2])
    lU <- rnorm(ng,mean=log(nU[1]),sd=nU[2])
    lTh <- lThU + lU
    lHe <- log(HeGain*IsoplotR:::get_He(tt,exp(lU),exp(lTh)))
    lUHe <- lU-lHe
    lThHe <- lTh-lHe
    vlUHe <- exp(-lU) + exp(-lHe)
    vlThHe <- exp(-lTh) + exp(-lHe)
    rho <- exp(-lHe)/(sqrt(exp(-lU)+exp(-lHe))*sqrt(exp(-lThHe)+exp(-lHe)))
    suv <- rho*sqrt(vlUHe*vlThHe)
    covmat <- matrix(0,nrow=2*ng,ncol=2*ng)
    diag(covmat) <- c(vlUHe,vlThHe)
    i1 <- 1:ng
    i2 <- (ng+1):(2*ng)
    i3 <- (2*ng+1):(3*ng)
    diag(covmat[i1,i2]) <- diag(covmat[i2,i1]) <- suv
    uv <- MASS::mvrnorm(1,mu=c(lUHe,lThHe),Sigma=covmat)
    u <- uv[i1]
    v <- uv[i2]
    den <- exp(u)+exp(v)+1
    J <- matrix(0,nrow=3*ng,ncol=2*ng)
    diag(J[i1,i1]) <- (exp(u)*den-exp(2*u))/den^2
    diag(J[i1,i2]) <- -exp(u+v)/den^2
    diag(J[i2,i1]) <- -exp(u+v)/den^2
    diag(J[i2,i2]) <- (exp(u)*den-exp(2*v))/den^2
    diag(J[i3,i1]) <- -exp(u)/den^2
    diag(J[i3,i2]) <- -exp(v)/den^2
    UThHe <- c(exp(u),exp(v),rep(1,ng))/den
    E <- J %*% covmat %*% t(J)
    out <- cbind(He=UThHe[i3]/HeGain,errHe=sqrt(diag(E)[i3])/HeGain,
                 U=UThHe[i1],errU=sqrt(diag(E)[i1]),
                 Th=UThHe[i2],errTh=sqrt(diag(E)[i2]),
                 Sm=0,errSm=0,t.true=tt)
    class(out) <- append(class(out),"UThHe")
    out
}
# converts a fission track age to log(rho_s/rho_i)
t2mu <- function(tt,zeta=350e-6,rhoD=2.5e6){
    l38 <- settings('lambda','U238')[1]
    log(2*(exp(l38*tt)-1)/(l38*zeta*rhoD))
}
# estimates approximate age uncertainty from age
t2cv <- function(tt,method='fissiontracks',
                 N=100,zeta=350e-6,rhoD=2.5e6,
                 ThU=1,nU=10000){
    if (identical(method,'fissiontracks')){
        mu <- t2mu(tt=tt,zeta=zeta,rhoD=rhoD)
        Ni <- N/(1+exp(mu))
        Ns <- N-Ni
        out <- sqrt(1/Ns+1/Ni)
    } else if (identical(method,'U-Th-He')){
        U <- nU
        Th <- ThU*U
        He <- IsoplotR:::get_He(tt,U,Th)
        sU <- sqrt(U)
        sTh <- sqrt(Th)
        sHe <- sqrt(He)
        tst <- IsoplotR::age(c(U,sU,Th,sTh,He,sHe),method='U-Th-He')
        out <- tst[2]/tst[1]
    } else {
        stop('Unrecognised method.')
    }
    out
}

dtrunc <- function(mux,mean,sd,minx=-Inf){
    deficit <- pnorm(minx,mean=mean,sd=sd)
    out <- dnorm(mux,mean=mean,sd=sd)/(1-deficit)
    out[mux<minx] <- 0
    out
}
rtrunc <- function(ng,mean,sd,minx=-Inf){
    out <- rnorm(ng,mean=mean,sd=sd)
    out[out<minx] <- minx[out<minx]
    out
}

# write fission track data to an IsoplotR input file
writeFTfile <- function(FTdata,fname){
    sink(fname)
    cat('zeta,err\n')
    cat(FTdata$zeta[1],',',FTdata$zeta[2],'\n',sep='')
    cat('rhoD,err\n')
    cat(FTdata$rhoD[1],',',FTdata$rhoD[2],'\n',sep='')
    cat('Ns,Ni\n')
    for (i in 1:length(FTdata)){
        cat(FTdata$x[i,1],',',FTdata$x[i,2],'\n',sep='')
    }
    sink()
}

# write U-Th-He data to an IsoplotR input file
writeUThHefile <- function(UThHedata,fname){
    write.csv(UThHedata[,1:6],fname,row.names=FALSE)
}
