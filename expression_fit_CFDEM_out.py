# -*- coding: utf-8 -*-
"""
Script to output the PSD for a Rosin-Rammler distribution
inputs:
  dMean: mean particle diameter
  dMin:  min particle diameter
  dMax:  max particle diameter
  n:  distribution coefficient
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import math

import json

from lmfit.models import ExpressionModel


#Rosin-Rammler Function taken from openfoam
#equation used from openfoam website for the fitting function...
#https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1distributionModels_1_1RosinRammler.html
def cfd_rr_OF(dMin,dMax,dMean,n,x):
  top = (1-np.exp((-x/dMean)**n+(dMin/dMean)**n))
  bottom = (1-np.exp(-(dMax/dMean)**n+(dMin/dMean)**n))
  return top/bottom

#Rosin-Rammler Function
def cdf_rr_func(x,a,b):
      return -np.power(x/a,b)

def rosin_rammler_CFDEM_out(fitted_variables, fitting_function,sample_name,dMin,dMax,num_particle_groups,particle_density,weight=1):

  particle_sizes = np.linspace(dMin,dMax,num_particle_groups)
  cdf = fit_function_out(particle_sizes,fitting_function,fitted_variables)

  #cdf_min is the reaming mass fraction below the minimum particle diameter
  cdf_min = cdf.min()
  print("the remaing mass fraction below dmin:" ,cdf_min)


  #scaling the cdf so that all mass fraction starts at 0 instead of dmin
  cdf_scaled = cdf - cdf_min
  cdf_scaled_max = cdf_scaled.max()

  print("cdf is:", cdf)

  print("modified cdf maximum:", cdf_scaled_max)

  cdf_scaled = cdf_scaled/cdf_scaled_max

  print("scaled cdf:", cdf_scaled)

  
  # print("the total cdf is: ", cdf)
  # print("the sume of the cdf is:", cdf.sum())

  #convert from a cummulative distribution to percentages of each particle size
  percentages = []
  last_index = cdf_scaled.size-1
  percentages.append(0)
  i = 1
  while i <= last_index:
    percentages.append(cdf_scaled[i]-cdf_scaled[i-1])
    i=i+1
    

  sum_of_percentages = np.sum(percentages)
  print("the percentages are:", percentages)
  print('sum of percentages:', sum_of_percentages)
  remaining_mass = 1 - sum_of_percentages

  percentages_cumsum = np.cumsum(percentages)

  print('percentages cumsum:', percentages_cumsum)

  file = open("CFDEM_output.txt","w+")

  i = 0
  for percentage in percentages:
    file.write(f'fix   pts{i} all particletemplate/sphere 1 atom_type 1 density constant {particle_density}\
    radius constant {particle_sizes[i]/2}e-6\n')
    i = i + 1

  file.write(f'fix pdd all particledistribution/discrete 63243 {num_particle_groups}')

  i = 0
  for percentage in percentages:
    file.write(f' pts{i} {percentage} ')
    i = i + 1

  file.close()

  xx_lagrangian = np.linspace(dMin,dMax,num=50)
  xx_eulerian = np.linspace(0,dMin,num=50)
  yvals_lagrangian=fit_function_out(xx_lagrangian,fitting_function,fitted_variables)
  yvals_eulerian=fit_function_out(xx_eulerian,fitting_function,fitted_variables)

  exp_data = np.genfromtxt(sample_name)
  exp_x = exp_data.transpose()[0]
  exp_y = exp_data.transpose()[1]

  #plot1=plt.plot(particle_sizes, percentages_cumsum, linewidth = 2,c='r',label='CFDEM_output')
  plot1=plt.plot(exp_x, exp_y, 'o',c='g',label='sample data')
  plot2=plt.plot(xx_lagrangian, yvals_lagrangian, 'b',linewidth=2,label='CFDEM_output (lagrangian)')
  plot3=plt.plot(xx_eulerian, yvals_eulerian, 'r',linewidth=2,label='eulerian fines')
  plt.xlabel('diameter(μm)')
  plt.ylabel('mass fraction')
  plt.legend(loc=1)
  plt.savefig('CFDEM_out.png')

  plt.xlabel('diameter(μm)')
  plt.ylabel('mass fraction')
  plt.legend(loc=1)

  return cdf_min


def fit_expression(sample_name,fitting_expression,**fitting_variables):
  data = np.genfromtxt(sample_name)
  x = data.transpose()[0]
  y = data.transpose()[1]
  x = x[:-1]
  y = y[:-1]

  gmod = ExpressionModel(fitting_expression)

  result = gmod.fit(y, x=x, **fitting_variables)

  print(result.fit_report())

  #plot5 = plt.plot(x, result.best_fit , 'p-',label=f'fit:{fitting_function}')

  return result


#using: https://www.cheresources.com/invision/blog/4/entry-340-calculating-physical-properties-of-slurries/
def eulerian_fines(remaining_mf,particle_mass_flow_rate,particle_density,fluid_mass_flow_rate,fluid_density,fluid_viscocity,dense=False):
  fines_particle_mass_flow_rate = particle_mass_flow_rate*remaining_mf

  print("lagrangian Mass Flow rate (minus fines): ",particle_mass_flow_rate-fines_particle_mass_flow_rate)
  print("fines (Eulerian) mass flow rate:", fines_particle_mass_flow_rate)

  mass_ratio = fines_particle_mass_flow_rate/fluid_mass_flow_rate
  mass_concentration = mass_ratio/(1+mass_ratio)
  print("mass ratio: ", mass_ratio)
  print("mass_concentration: ", mass_concentration)
  rho_m = 100/(mass_concentration*100/particle_density+((100-mass_concentration*100)/fluid_density))

  volume_concentration = mass_concentration*rho_m/particle_density

  print("volume concentration:", volume_concentration)

  if dense:
    mu_m = fluid_viscocity*(1 + 2.5*volume_concentration + 10.05*volume_concentration**2 + 0.00273*np.exp(16.6*volume_concentration))
  else:
    mu_m = fluid_viscocity*(1 + 2.5*volume_concentration)
  nu_m = mu_m/rho_m

  print("Fines Eulerian mixture density: ", rho_m)
  print("Fines Eulerian dynamic mixture viscocity: ", mu_m)
  print("Fines Eulerian kinematic mixture viscocity: ", nu_m)
  

def fit_function_out(x,fitting_expression,fitted_variables):

  math_dict = {'exp':np.exp, 'sin': np.sin}
  fitted_variables['x'] = x
 
  return eval(fitting_expression,fitted_variables,math_dict)


#input user inputs from .json file
with open('inputs.json') as f:
  data = json.load(f)

dMin = data["minimum_particle_dimater"]
dMax = data["maximum_particle_diameter"]
num_particle_groups = data["CFDEM_num_particle_groups"]
particle_density = data["particle_density"]
sample_name = data["sample_name"]
fitting_function = data["fitting_function"]
fitting_variables = data["fitting_variables"][0]

particle_mass_flow_rate = data["mass_flow_rate_particles"]
fluid_mass_flow_rate = data["mass_flow_rate_fluid"]
fluid_density = data["fluid_density"]
fluid_viscocity = data["fluid_dyn_viscocity"]

result = fit_expression(sample_name,fitting_function,**fitting_variables,method='leastsq')

fitted_variables = dict(result.params.valuesdict())


#plot the fitted equation with sample
x = np.linspace(0,dMax,num=50)
y = fit_function_out(x,fitting_function,fitted_variables)
exp_data = np.genfromtxt(sample_name)
exp_x = exp_data.transpose()[0]
exp_y = exp_data.transpose()[1]
plot1=plt.plot(exp_x, exp_y, 'o',c='g',label='sample data')
plt.plot(x, y, 'y-',label=f'CFDEM_output with fit:{fitting_function}')
plt.xlabel('Diameter(μm)')
plt.ylabel('Cumm% Passing')
plt.legend(loc=2)

plt.savefig('fitted_expression.png')

plt.close()


remaining_mf = rosin_rammler_CFDEM_out(fitted_variables, fitting_function,sample_name,dMin,dMax,num_particle_groups,particle_density)


eulerian_fines(remaining_mf,particle_mass_flow_rate,particle_density,fluid_mass_flow_rate,fluid_density,fluid_viscocity,dense=False)




