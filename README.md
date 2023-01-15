# Peak fitting of KIE measurement from Mass Spectrometry Data

This code has been applied on the following publications in order to accelerate the data processing efficiency of peak fitting and was verified by the traditional method with [Origin](https://www.originlab.com).

- **Hsuan-Chun Lin**, Benjamin Weissman, Syed Shahbaz Gardezi, Vernon Anderson, Darrin York, Joseph Piccirilli, Michael Harris **2018 ACS National meeting New Orleans, March 22-26**
Kinetic isotope effects on catalysis by the HDV ribozyme-precise determination of isotope ratios using electrospray ionization time-of-flight mass spectrometry

- Benjamin Weissman, Şölen Ekesan, **Hsuan-Chun Lin**, Shahbaz Gardezi, Nansheng Li, Timothy J. Giese, Erika McCarthy, Michael E Harris, Darrin M York, and Joseph A Piccirilli, A dissociative transition state in Hepatitis Delta Virus ribozyme catalysis. **J. Am. Chem. Soc.** in Press 2023.

![png](./images/Process.png)



```python
#import data
import pandas as pd
import numpy as np
import math
from scipy import asarray as ar
from scipy.optimize import leastsq
from scipy import integrate
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob, os
```

## Preprocessing data
The measurement of different concentration points are exported from machine and stored in .xy form.


```python
#Fetch file list
filelist = []
df = pd.DataFrame()
for file in glob.glob("*.xy"):
    filelist.append(file)
    print(file)
```

    RP2-0.xy
    RP2-90a.xy
    RP2-30.xy
    RP2-70.xy
    RP2-50.xy


## Define parameters (peaks)


```python
def generate_bounds(lb_a, lb_peak, lb_sigma, lb_baseline, peak_number):
    lb_as = [lb_a for i in range(peak_number)]
    lb_peaks = [lb_peak for i in range(peak_number)]
    lb_sigmas = [lb_sigma for i in range(peak_number)]
    exports = list(zip(lb_as,lb_peaks, lb_sigmas))
    exports = list(np.reshape(exports, -1))
    exports.insert(0, lb_baseline)
    return np.array(exports)
```


```python
#How many peaks you want to include
peak_number = 9

#Peak ranges on the spectrum
peak_range = (1156.0, 1159.0)

#Set initial values for fitting
#center of each peaks 
peak_center_s = [1156.1625, 1156.4968, 1156.8307, 1157.1649, 1157.4991, 1157.8332,
               1158.1668, 1158.5, 1158.8337]
#peak height -- use the same ones
a_s = [30000 for i in range(len(peak_center))]

#Gaussian sigmas -- also use identical ones
sigma_s = [0.1 for i in range(len(peak_center))]

#
baseline = 0

#Generate upper and lower bounds for baseline, as, peaks, peakhights
lower_bound = generate_bounds(0., 1155., 0., 0., 9)
upper_bound = generate_bounds(np.inf, 1159., 0.5, np.inf, 9)
f_bound = (lower_bound, upper_bound)
```

## Create Gaussian function for fitting
$$
baseline + ae^{- \frac{(x-x1)^2}{2 \sigma ^2}}
$$
baseline + a1*exp(-(x-x1)**2/(2*sigma1**2))

### Setting_fitting_function


```python
def string_for_fitting(peak_number):
    
    function_list = []
    function_vars = []
    for i in range(peak_number):
        component_a = '*np.exp(-(x-'
        component_b = ')**2/(2*'
        component_c = '**2))'
        V1='a' + str(i+1)
        V2='x' + str(i+1)
        V3='sigma'+str(i+1)
        string_gen = V1 +component_a + V2 + \
        component_b + V3 + component_c
        function_list.append(string_gen)
        #generate the function string
        export_str = ""
        #generate the function variable
        function_vars.append(V1)
        function_vars.append(V2)
        function_vars.append(V3)

    #generate function string
    for func in function_list:
            export_str += "+" + func   
    export_str = export_str[1:]
    
    #string for var
    start_string_var = "def gaus(x, baseline"
    for item in function_vars:
        start_string_var += "," + item
    start_string_var += "):\n\t"
    generate_function = start_string_var + "return(baseline + " + export_str + ")"

    return generate_function, export_str
```

## Fit the peaks


```python
def fit_peaks(x, y, peak_number):
    y_real = y
    import_data = list(zip(a_s, peak_center_s, sigma_s))
    import_data = list(np.reshape(import_data, -1))
    import_data.insert(0, baseline)

    popt,pcov = curve_fit(gaus,x,y,p0=import_data, bounds = f_bound)
    
    y_est = gaus(x, *popt)
    
    plt.plot(x, y_real,'g.',label='Real Data')
    plt.plot(x, y_est, label='Fitted')
    plt.legend()
    plt.show()

    return popt, pcov
```

## Integrate the peaks to calculate the area


```python
def integration(export_str, peak_number):
    inte_str= "inte = lambda x:" +  export_str
    exec(inte_str)
    #set variables
    for i in range(1,peak_number+1):
        comd1 = "x" + str(i) + "= peak"+str(i)
        comd2 = "a" + str(i) + "= a"+str(i)
        comd3 = "sigma" + str(i) + "= sigma"+str(i)

        print(comd1)
        print(comd2)
        print(comd3)

        exec(comd1)
        exec(comd2)
        exec(comd3)
    integration_result, integration_err = integrate.quad(inte, peak_range[0], peak_range[1] )
    print(integration_result, integration_err)
    return integration_result, integration_err
```

## Final calculation
### Generate Equation


```python
generate_function, export_str = string_for_fitting(peak_number)
exec(generate_function)
```

### Fitting loop


```python
for item in filelist:

    #For each .xy file, data preparation
    data = pd.read_table(item, sep = " ", header = None )
    fname = str(item)

    peaks = data[1].nlargest(5).values
    gauss1 = peaks[2]
    gauss2 = peaks[0]
    middle_point = data[data[1] == gauss1]
    middle_point2 = data[data[1] == gauss2]
    
    fittingData_1 = data[data[0] > peak_range[0]]
    fittingData_2 = fittingData_1[fittingData_1[0] < peak_range[1]]
    x = fittingData_2.iloc[:,0].values
    y = fittingData_2.iloc[:,1].values

    #Ready for fitting
    popt, _ = fit_peaks(x, y, peak_number)

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
        func = lambda x, a, sigma: a*np.exp(-(x)**2/(2*sigma**2))
        M1, M1_err = integrate.quad(func, -np.inf, np.inf, args = (a, sigma) )
        integration_array.append(M1)

    export_list = []
    for string in integration_array:
        export_list.append(string)

    df[fname] = export_list
```


    
![png](./images/output_16_0.png)
    



    
![png](./images/output_16_1.png)
    



    
![png](./images/output_16_2.png)
    



    
![png](./images/output_16_3.png)
    



    
![png](./images/output_16_4.png)
    


## Fitting Result
The result shows surface areas of peaks from 1 to 9 


```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>RP2-0.xy</th>
      <th>RP2-90a.xy</th>
      <th>RP2-30.xy</th>
      <th>RP2-70.xy</th>
      <th>RP2-50.xy</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>RP2-0.xy</td>
      <td>RP2-90a.xy</td>
      <td>RP2-30.xy</td>
      <td>RP2-70.xy</td>
      <td>RP2-50.xy</td>
    </tr>
    <tr>
      <th>1</th>
      <td>379.36016</td>
      <td>115.588081</td>
      <td>200.914163</td>
      <td>112.857305</td>
      <td>148.536372</td>
    </tr>
    <tr>
      <th>2</th>
      <td>518.374617</td>
      <td>159.810319</td>
      <td>274.409417</td>
      <td>153.840419</td>
      <td>202.418131</td>
    </tr>
    <tr>
      <th>3</th>
      <td>765.75728</td>
      <td>238.778048</td>
      <td>409.798355</td>
      <td>230.680275</td>
      <td>302.871698</td>
    </tr>
    <tr>
      <th>4</th>
      <td>707.76019</td>
      <td>220.810957</td>
      <td>379.401816</td>
      <td>213.308806</td>
      <td>279.14439</td>
    </tr>
    <tr>
      <th>5</th>
      <td>478.42302</td>
      <td>148.895244</td>
      <td>254.527133</td>
      <td>143.70668</td>
      <td>187.34727</td>
    </tr>
    <tr>
      <th>6</th>
      <td>253.283165</td>
      <td>80.35973</td>
      <td>137.074781</td>
      <td>78.03461</td>
      <td>101.893598</td>
    </tr>
    <tr>
      <th>7</th>
      <td>105.835663</td>
      <td>33.463338</td>
      <td>56.976236</td>
      <td>31.812005</td>
      <td>41.23961</td>
    </tr>
    <tr>
      <th>8</th>
      <td>34.310296</td>
      <td>11.079677</td>
      <td>18.26443</td>
      <td>10.23424</td>
      <td>13.176691</td>
    </tr>
    <tr>
      <th>9</th>
      <td>7.594068</td>
      <td>4.200387</td>
      <td>5.392131</td>
      <td>4.585807</td>
      <td>4.343414</td>
    </tr>
  </tbody>
</table>
</div>




```python
df.to_csv('fitting_result.csv')
```
