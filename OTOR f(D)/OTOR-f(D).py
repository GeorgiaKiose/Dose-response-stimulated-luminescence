#Fitting f(D) Vs Dose curves with the OTOR model
#Author: Georgia Kioselaki
#Contact details: gkioselaki@auth.gr

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.optimize as optimization
from scipy.special import lambertw
from tkinter import *
from tkinter import messagebox

'''Message for the user'''
window=Tk()
window.eval('tk::PlaceWindow %s center' % window.winfo_toplevel())
window.withdraw()
messagebox.showinfo("Message for the user", "Before you use the python script, remember to save the experimental data in a file with name 'data.txt' and the initial values for the parameters in a file 'params.txt'. Both files must have a txt (tab delimited) format. Keep them in the same folder with the python script!\n")
window.deiconify()
window.destroy()
window.quit()

'''Import of the experimental data and the initial parameter values that the user has defined.'''
with open ('data.txt', mode='r') as file:
	rawdata=list(csv.reader(file,delimiter='\t'))
data=np.array(rawdata[1:], dtype=float)
xdata=data[:,0]
ydata=data[:,1]

with open ('params.txt', mode='r') as file:
	rawdata2=list(csv.reader(file,delimiter='\t'))
i=0
guess=[]
for i in range(2):
     k = np.float (rawdata2 [i] [1])
     guess.append(k)
     i=i+1

'''Use of the analytical equation that describes the OTOR model.'''     
def fourPL2(x, R, Dc):
     ZR= (R-1)* np.exp (R-1)
     Z1= ZR *  np.exp (-x/Dc)
     W=np.real(lambertw(Z1))
     W2=np.real(lambertw(ZR))
     K= ((1/((R-1)*Dc)) * (W2/ (1+W2)))
     return (1/(K*x))*(1- ((W/(R-1))))

'''Fitting of the theoretical and the experimental curves. The parameter bounds are set.'''
params, params_covariance = optimization.curve_fit(fourPL2, xdata, ydata, guess, bounds= [0, ( 0.99999, np.inf)] ,maxfev=500000000)

#Output data
window=Tk()
window.eval('tk::PlaceWindow %s center' % window.winfo_toplevel())
window.withdraw()
messagebox.showinfo("Parameters", "The values of the fitting params are:\n R=%5.3f\n Dc=%5.3f.\n They have been saved in a .txt file in your folder." %tuple (params))

window.deiconify()
window.destroy()
window.quit()


'''FOM Calculation.'''
Ith=[]

for i in xdata:
	values=fourPL2(i, *params)
	Ith.append(values)
FOM = 100*(np.sum(abs(Ith-ydata))/(np.sum(ydata)))
print ('\nFOM=%0.3f' %FOM , '%')


'''The ascii file is written and includes the values of the theoretical curve'''
zip(xdata,Ith)
with open('Theoretical curve.txt', 'w', newline="") as f:
    writer = csv.writer(f,delimiter='\t')
    writer.writerow([rawdata [0] [0], "Ith (A.U.)"])
    writer.writerows (zip (xdata, Ith))


'''Creation of the txt file that includes the parameter values after the fitting.'''
parameters=['R=', 'Dc=']
with open('Fitting Params.txt', 'w', newline="") as f:
    writer= csv.writer(f, delimiter='\t')
    n=0
    for n in range(2):
         writer.writerow([parameters[n], params[n]])
         n=n+1

'''Graphical plotting and data visualization.'''
x_min, x_max = np.amin(xdata), np.amax(xdata)
xs = np.linspace(x_min, x_max, 100000)
plt.scatter(xdata, ydata, label='experimental data \nFOM=%5.3f' %FOM)

plt.plot(xs, fourPL2(xs, *params),color='y', label='fit params: \n R=%5.3f\n Dc=%5.3f\n' %tuple (params))

plt.legend()
plt.xlabel(rawdata [0] [0])
plt.ylabel(rawdata [0] [1])

plt.grid()
plt.xscale('log')
plt.savefig ('plot.png', bbox_inches='tight')
plt.show()


