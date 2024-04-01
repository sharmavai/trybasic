import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker
# Taking an semi-ellipse body for pressure distribution
#creating an semi-ellipse
t = np.linspace(0,np.pi,10000);
a = 1.5;
b = 0.2;
gamma = 1.4
#parametric equation of Ellipse
x_cor = 2 + a * np.cos(t); 
y_cor =  b * np.sin(t);
x_cord = np.flip(x_cor)
y_cord = np.flip(y_cor)
# Plot of the  body
plt.hlines(0, xmin=0, xmax=4)
plt.plot(x_cord,y_cord)
plt.xlim([0,4])
plt.ylim([-0.2,2])
plt.show()
theta_line = np.zeros(t.size-1)
for i in range(len(theta_line)):  #finding slope of the segments
  theta_line[i] = np.arctan( (y_cord[i+1]-y_cord[i])/(x_cord[i+1]-x_cord[i]))

#condering the flow is moving from left to right
M_inf = 6; #Taking mach of freee stream as 7
#For oblique shock theta beta mach relation
betamin = np.arcsin(1/M_inf)
beeta = np.array(np.arange(betamin,np.pi/2,np.pi/20000)) # storing values of beta for which theta would be calculated
theeta = np.arctan((2 * ((M_inf * np.sin(beeta))**2 - 1) / (np.tan(beeta)*(M_inf**2 * (1.4+np.cos(2*beeta)) + 2))))# using reltion of theta,beta and M
tmax =  np.max(theeta); #theta max 
max_position = pd.Series(theeta).idxmax() # finding index at which theta is maximun in theeata array
pressure = np.zeros(t.size-1); # 1 atm atmospheric pressure as flow static pressure
M = np.zeros(t.size-1)
#For expansion wave Prandtl–Meyer equation 
MACH_N = np.linspace(1,30,100000)
niu = np.sqrt((gamma+1)/(gamma-1)) * np.arctan(np.sqrt(((gamma-1)*(MACH_N**2-1))/(gamma+1)))-np.arctan(np.sqrt(MACH_N**2-1))
for i in range (0,theta_line.size-1):
  # creating a condition to check if the flow is subsonic or super sonic 
  # assumption of the normal shock throughout that subsonic portion untill theta approach theta max 
  if((theta_line[i]) >= tmax): 
    # another condition to conform that the speed changes from subsonic to super sonic after theta max
    if(np.absolute(theta_line[i]-tmax) < 0.0005235):
      beta1 =  beeta[0:int(max_position)]
      beta = beta1[-1]
      Mn = M_inf * np.sin(beta)
      pressure[i] =1 * (1+2*gamma*(Mn**2 - 1) / (gamma + 1)); 
      M1_n = np.sqrt((1+((gamma-1)/2)*Mn**2)/(gamma*Mn**2-(gamma-1)/2));
      M[i]=M1_n/(np.sin(beta-tmax))
    else:  
      Mn = M_inf ;
      pressure[i] = 1 * (1+2*gamma*(Mn**2 - 1) / (gamma + 1));   
      M1_n = np.sqrt((1+((gamma-1)/2)*Mn**2)/(gamma*Mn**2-(gamma-1)/2)); 
      M[i] =  M1_n;
      pressure[i] = pressure[i] * (1+((gamma-1)*M[i]**2)/2)**(gamma/(gamma-1)) ;
 #after sonic condition the expansion theory can be applied
  else:
    if(M[i-1]<1):
      M[i-1] = 1;
    nu = np.sqrt((gamma+1)/(gamma-1)) * np.arctan(np.sqrt(((gamma-1)*(M[i-1]**2-1))/(gamma+1)))-np.arctan(np.sqrt(M[i-1]**2-1))
    l_theta = theta_line[i-1]-theta_line[i]
    nu_2 = nu + l_theta
    difference_array = np.absolute(niu-nu_2)
    index = difference_array.argmin()
    M[i] = MACH_N[index] 
    pressure[i] = pressure[i-1] * ((1+((gamma-1)/2)*M[i-1]**2)/(1+((gamma-1)/2)*M[i]**2))**(gamma/(gamma-1))
# At the trailing edge slope suddenly increases hence the flow is turned into itself creating shock wave  
X_coordinate = np.zeros(t.size-1)
Y_coordinate = np.zeros(t.size-1)
# Assuming pressure is acting at the mid point of line segment
for i in range(0,len(x_cord)-1):
  X_coordinate[i] = ((x_cord[i+1]+x_cord[i])/2)
  Y_coordinate[i] = ((y_cord[i+1]+y_cord[i])/2)
cp = (pressure - 1)/(0.5*1.225*2324**2)*(101325)
plt.plot(X_coordinate,cp,linewidth = '1')
plt.xlabel('Chord (mm)')
plt.ylabel('Cp')
plt.title('Cp distribution')
plt.grid(color='b', linestyle='-', linewidth=0.1)
plt.xlim([0,4])
plt.show() 
#for cp location 
#xpdx =0;
#pdx = 0
#for i in range (0,x_cord.size-1):
# xpdx = xpdx + (X_coordinate[i]*pressure[i]*(x_cord[i+1]-x_cord[i])+ X_coordinate[i] * 1 * (x_cord[i+1]-x_cord[i]) *np.cos(theta_line[i]))*0.001*101325
#  pdx = pdx + pressure[i]*(x_cord[i+1]-x_cord[i])*101325*0.001

#Xcp = xpdx/pdx # cp location 





