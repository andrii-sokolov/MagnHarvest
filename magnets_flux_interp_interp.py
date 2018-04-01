import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return(a**2/(b**2+x**2)**1.5)

def func_a(r,c,d,e,f,g):
    return(c + d*r + e*r**2 + f*r**3 + g*r**4)

def FluxFile():
    file = open('/home/cbapka/cpp_flux_to_interp.dat','r')
    
    ## Transforming data from files to float numbers
    data1 = []
    for line in file:
        data1.append(line.split())
    data2 = []
    for d in data1:
        temp = []
        for d1 in d:
            temp.append(float(d1))
        data2.append(temp)

    data3 = []
    for i in range(len(data2[0])):
        data3.append([])
    for d in data2:
        for i in range(len(d)):
            data3[i].append(d[i])
    return(data3)

fitz_r=[]
fitz_a=[]
fitz_b=[]

dat = FluxFile()
for i in range(1,len(dat)):
    popt, pcov = curve_fit(func, dat[0], dat[i])
    fitz_r.append(0.01+0.001*i)
    fitz_a.append(np.fabs(popt[0]))
    fitz_b.append(np.fabs(popt[1]))

popt_1, pcov_1 = curve_fit(func_a, fitz_r, fitz_a)
popt_2, pcov_2 = curve_fit(func_a, fitz_r, fitz_b)

print(popt_1)
print(popt_2)

fited=[]
for r in fitz_r:
    fited.append(func_a(r,popt_1[0],popt_1[1],popt_1[2],popt_1[3],popt_1[4]))

plt.plot(fitz_r, fited)
plt.plot(fitz_r, fitz_a,'o')

plt.ylabel('a')
plt.xlabel('r, m')
plt.show()

fited=[]
for r in fitz_r:
    fited.append(func_a(r,popt_2[0],popt_2[1],popt_2[2],popt_2[3],popt_2[4]))

plt.plot(fitz_r, fited)
plt.plot(fitz_r, fitz_b,'o')

plt.ylabel('b')
plt.xlabel('r, m')
plt.show()


fitted=[]
for d in np.arange(-0.075,0.075,0.001):
    fitted.append(func(d,func_a(0.01+0.001*20,popt_1[0],popt_1[1],popt_1[2],popt_1[3],popt_1[4]),func_a(0.01+0.001*20,popt_2[0],popt_2[1],popt_2[2],popt_2[3],popt_2[4])))

plt.plot(dat[0], dat[20],'o')
plt.plot(np.arange(-0.075,0.075,0.001), fitted)

plt.ylabel('Phi, Wb')
plt.xlabel('r, m')
plt.show()
