############################################################
##a simple model for different states of a protein in a signaling cascade.
##written by Mengyue Lei (1155105728@link.cuhk.edu.hk), last change 1/15/2018
############################################################
rm(list=ls()) #this clears the workspace to make sure no leftover variables are floating around. Not strictly needed
graphics.off(); #close all graphics windows
library(deSolve)  #loads ODE solver package

#functions come first, main program below

###################################################################
#function that specificies the ode model called by lsoda (the ode solver) 
###################################################################
odeequations=function(t,y,pars) 
{ 
  #Note: y is a vector containing the variables of our system, pars is a vector containing the parameters
  P0=y[1]; P1=y[2]; P2=y[3]  #different states of a protein
  r01=pars[1]; r12=pars[2]; r21=pars[3]; r10=pars[4]; #model parameters, passed as vector "par" into function by main program
  
  #these are the differential equations
  dP0dt=r10*P1-r01*P0;
  dP1dt=r01*P0+r21*P2-r12*P1-r10*P1;
  dP2dt=r12*P1-r21*P2;
  
  return(list(c(dP0dt,dP1dt,dP2dt))); #this is returned to the calling function.
  
} #end function specifying the ODEs

###################################################################
#main program
###################################################################
#initial number
P00=2;  
P10=5;
P20=3;
Y0=c(P00, P10, P20);  #combine initial conditions into a vector 

#values for model parameters, units are assumed to be 1/hours
r01=0.3; 
r12=0.6;
r21=0.6;
r10=0.3;
pars=c(r01,r12,r21,r10); #vector of parameters which is sent to the ODE function

tmax=10; #number of days for which to run the simulation
timevec=seq(0,tmax,0.1); #vector of times for which integration is evaluated (from 0 to 10 hours in steps of 0.01)

#call ode-solver to integrate ODEs
#integrate for time "timevec", starting with initial condition 'Y0'. 
odeoutput=lsoda(Y0,timevec,odeequations,pars);

#plot results
#first column contains time vector, the following columns contain variables 1 (bacteria) and 2 (immune response)
plot(odeoutput[,1],odeoutput[,2],type="l",xlab="time(hours)",ylab="",col="blue",lwd=2, xlim=c(0,tmax),ylim=c(1,6))
lines(odeoutput[,1],odeoutput[,3],type="l",col="red",lwd=2)
lines(odeoutput[,1],odeoutput[,4],type="l",col="yellow",lwd=2)
legend("topright", c("P0","P1","P2"),col = c("blue","red","yellow"),lwd=2, cex = 0.5)

###################################################################
#end main program
################################################################### 
