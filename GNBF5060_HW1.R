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
  a=y[1]; b=y[2]; c=y[3]; d=y[4]  #different states of a protein
  k1=pars[1]; k2=pars[2]; k3=pars[3]; k4=pars[4]; k5=pars[5]; k6=pars[6];#model parameters, passed as vector "par" into function by main program
  
  #these are the differential equations
  dadt=k1-k2*a;
  dbdt=k2*a-k3*b;
  dcdt=k2*a+k5*d^2-k4*c;
  dddt=2*k4*c-k5*d^2-k6*d;
  
  
  return(list(c(dadt,dbdt,dcdt,dddt))); #this is returned to the calling function.
  
} #end function specifying the ODEs

###################################################################
#main program
###################################################################
#initial number
a=0;  
b=0;
c=0;
d=0;
Y0=c(a,b,c,d);  #combine initial conditions into a vector 

#values for model parameters, units are assumed to be 1/hours
k1=2; 
k2=1;
k3=3;
k4=2;
k5=2;
k6=6;
pars=c(k1,k2,k3,k4,k5,k6); #vector of parameters which is sent to the ODE function

tmax=100; #number of days for which to run the simulation
timevec=seq(0,tmax,0.1); #vector of times for which integration is evaluated (from 0 to 10 hours in steps of 0.01)

#call ode-solver to integrate ODEs
#integrate for time "timevec", starting with initial condition 'Y0'. 
odeoutput=lsoda(Y0,timevec,odeequations,pars);

#plot results
#first column contains time vector, the following columns contain variables 1 (bacteria) and 2 (immune response)
plot(odeoutput[,1],odeoutput[,2],type="l",xlab="time(hours)",ylab="",col="blue",lwd=2, xlim=c(0,tmax),ylim=c(0,5))
lines(odeoutput[,1],odeoutput[,3],type="l",col="red",lwd=2)
lines(odeoutput[,1],odeoutput[,4],type="l",col="yellow",lwd=2)
lines(odeoutput[,1],odeoutput[,5],type="l",col="black",lwd=2)
legend("topright", c("a","b","c","d"),col = c("blue","red","yellow","black"),lwd=2, cex = 0.5)

###################################################################
#end main program
################################################################### 