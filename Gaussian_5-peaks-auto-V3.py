# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 02:31:27 2017

@author: fatmi
"""

#import data
import pandas as pd
import numpy as np
import math
from scipy import asarray as ar,exp
from scipy.optimize import leastsq
from scipy import integrate
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob, os
#Process data
filelist = []
df = pd.DataFrame()
for file in glob.glob("*.xy"):
    filelist.append(file)
    print(file)
    # get the file list
    
    #start to analysis
for item in filelist:
    data = pd.read_table(item, sep = " ", header = None )
    fname = str(item)
    peak_number = 9
    peak_range = (1156.0, 1159.0)
    peaks = data[1].nlargest(5).values
    gauss1 = peaks[2]
    gauss2 = peaks[0]
    
    middle_point = data[data[1] == gauss1]
    middle_point2 = data[data[1] == gauss2]
    
    fittingData_1 = data[data[0] > peak_range[0]]
    fittingData_2 = fittingData_1[fittingData_1[0] < peak_range[1]]
    x = fittingData_2.iloc[:,0].values
    y = fittingData_2.iloc[:,1].values
    
    peak1 = 1156.1625
    peak2 = 1156.4968
    peak3 = 1156.8307
    peak4 = 1157.1649
    peak5 = 1157.4991
    peak6 = 1157.8332
    peak7 = 1158.1668
    peak8 = 1158.5
    peak9 = 1158.8337
    a1 = 30000
    a2 = 30000
    a3 = 30000
    a4 = 30000
    a5 = 30000
    a6 = 30000
    a7 = 30000
    a8 = 30000
    a9 = 30000
    sigma1 = 0.1
    sigma2 = 0.1
    sigma3 = 0.1
    sigma4 = 0.1
    sigma5 = 0.1
    sigma6 = 0.1
    sigma7 = 0.1
    sigma8 = 0.1
    sigma9 = 0.1
    baseline = 0
    #define bound
    lower_bound = (0,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.,0.,1155.,0.)
    upper_bound = (np.inf,np.inf,1159.,0.5,np.inf,1159.,0.5,np.inf,1159.,0.5,np.inf,1159.,0.5,np.inf,1159.,0.5,np.inf,1159.,0.5,np.inf,1159.,0.5,np.inf,1159.,0.3,np.inf,1159.,0.3)
    f_bound = (lower_bound, upper_bound)
    ######################################
    
    #setting_fitting_function
    myfunction = ""
    function_list = []
    function_vars = []
    for i in range(1,peak_number+1):
        component_a = '*exp(-(x-'
        component_b = ')**2/(2*'
        component_c = '**2))'
        V1='a' + str(i)
        V2='x' + str(i)
        V3='sigma'+str(i)
        string_gen = V1 +component_a + V2 + \
        component_b + V3 + component_c
        function_list.append(string_gen)
        #generate the function string
        export_str = ""
        #generate the function variable
        function_vars.append(V1)
        function_vars.append(V2)
        function_vars.append(V3)
    #string for function
    for func in function_list:
            export_str += "+" + func   
    export_str = export_str[1:]
    
    #string for var
    start_string_var = "def gaus(x, baseline"
    for item in function_vars:
        start_string_var += "," + item
    start_string_var += "):\n\t"
    generate_function = start_string_var + "return(baseline + " + export_str + ")"
    exec(generate_function)
    ######################################
    # Solving
    y_real = y
    
    #prepare import_data
    string_list = "import_data=[baseline,"
    for i in range(1,peak_number+1):
        string_list += 'a'+str(i) + ',peak' + str(i) + ",sigma" + str(i) + ","
        
    string_list = string_list[:-1] + "]"
    print(string_list)
    exec(string_list)   
    #import_data = [a1,peak1,sigma1,a2,peak2,sigma2,a3,peak3,sigma3,a4,peak4,sigma4,a5,peak5,sigma5]
    popt,pcov = curve_fit(gaus,x,y,p0=import_data, bounds = f_bound)
    
    y_est = gaus(x, *popt)
    
    plt.plot(x, y_real,'g.',label='Real Data')
    plt.plot(x, y_est, label='Fitted')
    plt.legend()
    plt.show()
    
    #integration
    inte_str= "inte = lambda x:" +  export_str
    exec(inte_str)
    #set variables
    for i in range(1,peak_number+1):
        comd1 = "x" + str(i) + "= peak"+str(i)
        comd2 = "a" + str(i) + "= a"+str(i)
        comd3 = "sigma" + str(i) + "= sigma"+str(i)
        exec(comd1)
        exec(comd2)
        exec(comd3)
    integration_result, integration_err = integrate.quad(inte, peak_range[0], peak_range[1] )
    #print(integration_result, integration_err)
    
    parameters = []
    for i in range(0,peak_number):
        pp = []
        for j in range(0,3):
            pp.append(popt[j+i*3+1])
        parameters.append(pp)
    
    #set single gaussian
    integration_array = [item]
    for item in parameters:
        a = item[0]
        sigma = item[2]
        x = item[1]
        func = lambda x, a, sigma: a*exp(-(x)**2/(2*sigma**2))
        M1, M1_err = integrate.quad(func, -np.inf, np.inf, args = (a, sigma) )
        integration_array.append(M1)
    export_list = []
    for string in integration_array:
        export_list.append(string)
    df[fname] = export_list
    #export to CSV
df.to_csv('fitting_result.csv')