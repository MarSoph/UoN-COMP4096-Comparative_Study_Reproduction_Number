library(DEoptim)
library(Metrics)
library(dplyr)
library(deSolve)
library(ggplot2)

SDpar=read.csv('covid-stringency-index.csv')
SDpar=SDpar %>% filter(Entity=='Cyprus', Day>='2021-01-01', Day<=	'2021-05-31')

Optimisation=function(v_value, niter, strat){
  SEPIJARV_opt=function(t,y,parms){
    a1_inv=parms[1]
    a1=1/a1_inv
    a2_inv=parms[2]
    a2=1/a2_inv
    k_inv=parms[3]
    k=1/k_inv
    kv=k
    vaccdelay_inv=parms[4]
    vaccdelay=1/vaccdelay_inv
    d=parms[5]
    n_inv=parms[6]
    n=1/n_inv
    nv=n
    f=parms[7]
    p=parms[8]
    pv=parms[9]
    u=parms[10]
    uv=u
    uA=parms[11]
    g=parms[12]
    gv=g
    l_inv=parms[13]
    l=1/l_inv
    lv=l
    v=v_value
    b=parms[14]
    b1=b
    db2=parms[15]
    b2=db2*b1
    db3=parms[16]
    b3=db3*b1
    db4=parms[17]
    b4=db4*b2
    db5=parms[18] 
    b5=db5*b1
    db6=parms[19]
    b6=db6*b2
    ro=(SDpar$stringency_index)/100
    with(as.list(y),{
      
      dS=-(1-ro[t+1])*b*S*(I+d*A)-(1-ro[t+1])*b2*S*(Iv+d*Av)-S*(b5*J+b6*Jv+b3*P+b4*Pv+v)
      dE=(1-ro[t+1])*b*S*(I+d*A)+(1-ro[t+1])*b2*S*(Iv+d*Av)+S*(b5*J+b6*Jv+b3*P+b4*Pv)+
        (1-ro[t+1])*b*Tv*(I+d*A)+(1-ro[t+1])*b2*Tv*(Iv+d*Av)+Tv*(b5*J+b6*Jv+b3*P+b4*Pv)-l*E
      dP=l*E-k*P
      dI=p*k*P-(1-u)*a1*I-u*g*I
      dJ=u*g*I-a2*J+uA*g*A
      dA=(1-p)*k*P-(1-uA)*n*A-uA*g*A
      
      dT=S*v-(1-ro[t+1])*b*Tv*(I+d*A)-
        (1-ro[t+1])*b2*Tv*(Iv+d*Av)- Tv*(b5*J+b6*Jv+b3*P+b4*Pv)-vaccdelay*Tv
      dV=vaccdelay*Tv-(1-ro[t+1])*b*V*(I+d*A)-(1-ro[t+1])*b2*V*(Iv+d*Av)-V*(b5*J+b6*Jv+b3*P+b4*Pv)
      dEv=(1-ro[t+1])*b*V*(I+d*A)+(1-ro[t+1])*b2*V*(Iv+d*Av)+V*(b5*J+b6*Jv+b3*P+b4*Pv)-lv*Ev
      dPv=lv*Ev-kv*Pv
      dIv=pv*kv*Pv-(1-uv)*a1*Iv-uv*gv*Iv
      dJv=uv*gv*Iv-a2*Jv
      dAv=(1-pv)*kv*Pv-n*Av
      dR=(1-uA)*n*A+(1-u)*a1*I+f*a2*J + n*Av+(1-uv)*a1*Iv+a2*Jv
      dD=(1-f)*a2*J
      
      Dout=c(dS,dE,dP,dI,dJ,dA,dT,dV,dEv,dPv,dIv,dJv,dAv,dR,dD)
      #print(sum(Dout))
      return(list(Dout))
    })
  }
  
  ode_MSE=function(x, y0=y0, data=Simu){
    a1_inv=x[1]
    a1=1/a1_inv
    a2_inv=x[2]
    a2=1/a2_inv
    k_inv=x[3]
    k=1/k_inv
    vaccdelay_inv=x[4]
    vaccdelay=1/vaccdelay_inv
    d=x[5]
    n_inv=x[6]
    n=1/n_inv
    f=x[7]
    p=x[8]
    pv=x[9]
    u=x[10]
    uv=u
    uA=x[11]
    g=x[12]
    l_inv=x[13]
    l=1/l_inv
    b=x[14]
    b1=b
    db2=x[15]
    b2=db2*b1
    db3=x[16]
    b3=db3*b1
    db4=x[17]
    b4=db4*b2
    db5=x[18] 
    b5=db5*b1
    db6=x[19]
    b6=db6*b2
    v=v_value
    
    if (v_value==0.013){
      Simu=read.csv('SimuNoisy0_013.csv')
    }
    if (v_value==0.005){
      Simu=read.csv('SimuNoisy0_005.csv')
    }
    if (v_value==0.025){
      Simu=read.csv('SimuNoisy0_025.csv')
    }
    
    E=7973/888000 #01/01/2021
    I=280/888000
    R=6845/888000
    Tv=0
    A=1200/888000
    V=0
    Ev=0
    Iv=0
    Av=0
    Jv=0
    J=2520/888000
    P=4000/888000
    Pv=0
    D=127/888000
    N=888000/888000
    y0=c(S=N-E-P-I-J-A-Tv-V-Ev-Pv-Iv-Jv-Av-R-D,E=E,
         P=P, I=I,J=J, A=A, Tv=Tv, V=V, Ev=Ev,
         Pv=Pv, Iv=Iv, Jv=Jv, Av=Av, R=R, D=D)
    
    out=ode(y=y0,times=seq(0,150, by=1),func = SEPIJARV_opt, parms=x, method='ode45')
    out=as.data.frame(out)
    rep=numeric(151)
    for (i in 1:151){
      rep[i+1]=u*g*out$I[i]+uA*g*out$A[i]+uv*g*out$Iv[i]
    }
    error=rmse(Simu$x, rep)
    return(error)
  }

  lower=c(3,5,1,7,0.3,3,0.8,0,0,0,0,0,1,0.2,0,0,0,0,0)
  upper=c(7,15,7,15,1,9,1,1,1,1,1,1,3,5,1,1,1,1,1)
  opt=DEoptim(ode_MSE ,lower,upper, control = list(strategy =strat ,itermax=niter,NP=200))
  
  

  optimal_parms=opt[['optim']][["bestmem"]]
  names(optimal_parms)=c('a1_inv','a2_inv','k_inv','vaccdelay_inv','d','n_inv','f','p','pv','u','uA','g','l_inv','b','b2','b3','b4','b5','b6')
  new_dat=as.data.frame(ode(y=y0,times=seq(0,150, by=1),
                            func = SEPIJARV_opt, 
                            parms=optimal_parms, 
                            method='ode45'))

  return(list(opt,new_dat,optimal_parms)) 
}

DE_calc_Rt=function(v_value,data, x){
  R0=numeric(151)
  for (t in 1:151){
    a1_inv=x[1]
    a1=1/a1_inv
    a2_inv=x[2]
    a2=1/a2_inv
    k_inv=x[3]
    k=1/k_inv
    kv=k
    vaccdelay_inv=x[4]
    vaccdelay=1/vaccdelay_inv
    d=x[5]
    n_inv=x[6]
    n=1/n_inv
    f=x[7]
    p=x[8]
    pv=x[9]
    u=x[10]
    uv=u
    uA=x[11]
    g=x[12]
    gv=g
    l_inv=x[13]
    l=1/l_inv
    lv=l
    b=x[14]
    b1=b
    db2=x[15]
    b2=db2*b1
    db3=x[16]
    b3=db3*b1
    db4=x[17]
    b4=db4*b2
    db5=x[18] 
    b5=db5*b1
    db6=x[19]
    b6=db6*b2
    v=v_value
    ro=(SDpar$stringency_index)/100
    
    parameters=list(S=1,E=0,P=0, I=0,J=0, A=0, 
                           Tv=0, V=0, Ev=0, Pv=0, 
                           Iv=0, Jv=0, Av=0, R=0, D=0, 
                       a1=a1, a2=a2, k=k,kv=kv,
                       d=d, n=n,  f=f, p=p, pv=pv,
                       vaccdelay=vaccdelay, u=u, uv=uv,
                       uA=uA, g=g, gv=gv, 
                       l=l, lv=lv, v=v_value,
                       ro=ro[t], b1=b, b2=b2, b3=b3, b4=b4,
                       b5=b5, b6=b6)
    
    f=with(parameters, matrix(c(eval(fEE),eval(fEP),eval(fEI),eval(fEJ),eval(fEA), 
                                eval(fEEv),eval(fEPv),eval(fEIv),eval(fEJv),eval(fEAv), 
                                eval(fPE),eval(fPP),eval(fPI),eval(fPJ),eval(fPA),
                                eval(fPEv),eval(fPPv),eval(fPIv),eval(fPJv),eval(fPAv),
                                eval(fIE),eval(fIP),eval(fII),eval(fIJ),eval(fIA),
                                eval(fIEv),eval(fIPv),eval(fIIv),eval(fIJv),eval(fIAv),
                                eval(fJE),eval(fJP),eval(fJI),eval(fJJ),eval(fJA),
                                eval(fJEv),eval(fJPv),eval(fJIv),eval(fJJv),eval(fJAv),
                                eval(fAE), eval(fAP),eval(fAI),eval(fAJ),eval(fAA),
                                eval(fAEv), eval(fAPv),eval(fAIv),eval(fAJv),eval(fAAv),
                                eval(fEvE),eval(fEvP),eval(fEvI),eval(fEvJ),eval(fEvA), 
                                eval(fEvEv),eval(fEvPv),eval(fEvIv),eval(fEvJv),eval(fEvAv), 
                                eval(fPvE),eval(fPvP),eval(fPvI), eval(fPvJ),eval(fPvA),
                                eval(fPvEv),eval(fPvPv),eval(fPvIv),eval(fPvJv),eval(fPvAv),
                                eval(fIvE),eval(fIvP),eval(fIvI),eval(fIvJ),eval(fIvA),
                                eval(fIvEv),eval(fIvPv),eval(fIvIv),eval(fIvJv),eval(fIvAv),
                                eval(fJvE), eval(fJvP),eval(fJvI),eval(fJvJ),eval(fJvA),
                                eval(fJvEv),eval(fJvPv),eval(fJvIv),eval(fJvJv),eval(fJvAv),
                                eval(fAvE), eval(fAvP),eval(fAvI),eval(fAvJ), eval(fAvA),
                                eval(fAvEv), eval(fAvPv),eval(fAvIv),eval(fAvJv),eval(fAvAv)),
                              nrow=10, byrow=TRUE))
    
    vv=with(parameters, matrix(c(eval(vEE),eval(vEP),eval(vEI),eval(vEJ),eval(vEA),
                                 eval(vEEv),eval(vEPv),eval(vEIv),eval(vEJv),eval(vEAv),
                                 eval(vPE),eval(vPP),eval(vPI),eval(vPJ),eval(vPA),
                                 eval(vPEv),eval(vPPv),eval(vPIv),eval(vPJv),eval(vPAv),
                                 eval(vIE),eval(vIP),eval(vII),eval(vIJ),eval(vIA),
                                 eval(vIEv),eval(vIPv),eval(vIIv),eval(vIJv),eval(vIAv),
                                 eval(vJE),eval(vJP),eval(vJI),eval(vJJ),eval(vJA),
                                 eval(vJEv),eval(vJPv),eval(vJIv),eval(vJJv),eval(vJAv),
                                 eval(vAE), eval(vAP),eval(vAI),eval(vAJ),eval(vAA),
                                 eval(vAEv), eval(vAPv),eval(vAIv),eval(vAJv),eval(vAAv),
                                 eval(vEvE),eval(vEvP),eval(vEvI),eval(vEvJ),eval(vEvA),
                                 eval(vEvEv),eval(vEvPv),eval(vEvIv),eval(vEvJv),eval(vEvAv),
                                 eval(vPvE),eval(vPvP),eval(vPvI), eval(vPvJ),eval(vPvA),
                                 eval(vPvEv),eval(vPvPv),eval(vPvIv),eval(vPvJv),eval(vPvAv),
                                 eval(vIvE),eval(vIvP),eval(vIvI),eval(vIvJ),eval(vIvA),
                                 eval(vIvEv),eval(vIvPv),eval(vIvIv),eval(vIvJv),eval(vIvAv),
                                 eval(vJvE),eval(vJvP),eval(vJvI),eval(vJvJ),eval(vJvA),
                                 eval(vJvEv),eval(vJvPv),eval(vJvIv),eval(vJvJv),eval(vJvAv),
                                 eval(vAvE), eval(vAvP),eval(vAvI),eval(vAvJ),eval(vAvA),
                                 eval(vAvEv),eval(vAvPv),eval(vAvIv),eval(vAvJv),eval(vAvAv)),
                               nrow=10, byrow=TRUE))
    
    
    R0[t]=(max(eigen(f %*% solve(vv))$values))
  }
  
  Rt=(data$S + data$Tv+data$V)*R0
  Rt=as.data.frame(Rt)
  Rt$t=seq(1:151)
  return(list(R0,Rt))
}


