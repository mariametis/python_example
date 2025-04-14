import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import scipy
from scipy import optimize as op
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import axes3d
from numpy import arange

from scipy.interpolate import splprep
#https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
#https://docs.scipy.org/doc/scipy/reference/interpolate.html
#https://stackoverflow.com/questions/18962175/spline-interpolation-coefficients-of-a-line-curve-in-3d-space
#https://stackoverflow.com/questions/42658190/plotting-one-sigma-error-bars-on-a-curve-fit-line-in-scipy
#https://datagy.io/python-return-multiple-values/

# define the true objective function
def objective(x, a, b, c, f):
 return (a * x) + (b * x**2) + (c * x**3)  + f

def residuals(x, a, b, c, f, y):     #for least square
      e=x*0+1.
      return (y-objective(x, a, b, c, f)) / e

def chisq(x, a, b, c,  f,y):
    return (residuals(x, a, b, c, f,y)**2).sum()

# define the true objective function
def objective4(x, a, b, c,d, f):
 return (a * x) + (b * x**2) + (c * x**3)  + + (c * x**4)+f

def residuals4(x, a, b, c,d, f, y):     #for least square
    e=x*0+1.
    return (y-objective4(x, a, b, c,d, f)) / e

def chisq4(x, a, b, c, d, f,y):
    return (residuals4(x, a, b, c,d, f,y)**2).sum()

def get_fit(xi2,zi2,x_line):
  popt, pcov = curve_fit(objective, xi2, zi2) 
  a, b, c,  f = popt
  y_line = objective(x_line, a, b, c,  f)
  yf=objective(xi2, a, b, c,  f)
  #e=xi2*0+1.
  resid=residuals(xi2, a, b, c,  f,zi2)
  stand=resid.std()
  #aa=yf-zi2
  #stand2=aa.std()
  #chi2=chisq(xi2, a, b, c,  f,zi2)
  print(popt)
  print(a,b,c,f)
  print(np.sqrt(np.diag(pcov)))
  return { 'Yfit': yf, 'y_line': y_line, 'std':stand}

def get_fit4(xi2,zi2,x_line):
  popt, pcov = curve_fit(objective4, xi2, zi2) 
  a, b, c, d, f = popt
  y_line = objective4(x_line, a, b, c,d,  f)
  yf=objective4(xi2, a, b, c,d,  f)
  resid=residuals(xi2, a, b, c,  f,zi2)
  stand=resid.std()
  return { 'Yfit': yf, 'y_line': y_line, 'std':stand}

def get_plot(xi2,zi2,yf,x_line, y_line,xtitle,ytitle):
  fig = plt.figure() 
  ax = fig.add_subplot(111) 
  ax.scatter(xi2, yf, color='cyan') 
  ax.scatter(xi2, zi2, color='pink') 
  ax.set_xlabel(xtitle) 
  ax.set_ylabel(ytitle) 
  ax.scatter(x_line, y_line, color='r') 
  plt.show()


#data = pd.read_csv('grid_rowMSX.tab', header=0 ,sep='\s+')
data = pd.read_csv('data_WISE_PLUS.csv')
print(data)

#######2nd iteration

data2=data[ (data['KW2der'].values >-100)& (data['KW3der'].values >-100)& (data['KW4der'].values >-100)] #- yf.values < 1.5) & (data['JKder'].values > 0 )]

############################
xi2=data2['qw3']
yi2=data2['qw3']
zi2=data2['KW3der']

x_line = arange(-25., 0., 0.01)
items=get_fit(xi2,zi2,x_line)
yf=items.get("Yfit")
y_line=items.get("y_line")
sigma=items.get("std")
KW3q=y_line
print(sigma)
get_plot(xi2,zi2,yf,x_line, y_line,'Qw3','(K-W3)o')
##############################################################

############################
xi2=data2['qw4']
yi2=data2['qw4']
zi2=data2['KW4der']

items=get_fit(xi2,zi2,x_line)
yf=items.get("Yfit")
y_line=items.get("y_line")
sigma=items.get("std")
qw34q=np.array(x_line)
KW4q=y_line
print(sigma)
get_plot(xi2,zi2,yf,x_line, y_line,'Qw4','(K-W4)o')

##############################################################
grid2=pd.DataFrame(pd.Series(qw34q))
grid2['qw34q']=pd.Series(qw34q)
grid2['KW3q']=pd.Series(KW3q)
grid2['KW4q']=pd.Series(KW4q)
grid2.to_csv('grid2_WISE_PLUSQ.csv')

##############################################################

par = pd.read_csv('par21Fritz.tab')
cW3den=(1.-par['aw3ak'].values)
XX=data2['qw3'].values
KW3fit=-0.1354624650055342*XX+0.0025538763327201614*XX**2+4.363383019080156e-05*XX**3+0.21997033514232298
Ak_w3=(data2['KW3obs']-KW3fit)/cW3den

cW4den=(1.-par['aw4ak'].values)
XX=data2['qw4'].values
KW4fit=-0.08871312859744133*XX+0.008130561313315675*XX**2+0.00013519614614686287*XX**3+0.7126895654171297
Ak_w4=(data2['KW4obs'].values-KW4fit)/cW4den

data2['KW3fit']=KW3fit
data2['Ak_w3']=Ak_w3
data2['KW4fit']=KW4fit
data2['Ak_w4']=Ak_w4
data2.to_csv('data_WISE_PLUSout.csv')



