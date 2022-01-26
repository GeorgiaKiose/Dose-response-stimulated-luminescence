#Fitting TL/OSL Dose response curves with the TTOR model
#Author: Georgia Kioselaki
#Contact details: gkioselaki@auth.gr

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.optimize as optimization
import csv
from scipy.special import lambertw
from tkinter import *
from tkinter import messagebox

'''Message for the user'''
window=Tk()
window.eval('tk::PlaceWindow %s center' % window.winfo_toplevel())
window.withdraw()
messagebox.showinfo("Message for the user", "Before you use the python script, remember to save the experimental data in a file with name 'data' and the initial values for the parameters in a file 'params'. Both files must have a text (tab delimited) format. Keep them in the same folder with the python script!")

window.deiconify()
window.destroy()
window.quit()
      

'''Import of the experimental data and the initial parameter values that the user has defined.'''
with open ('data.txt', mode='r') as file:
	rawdata=list(csv.reader(file,delimiter='\t'))
data=np.array(rawdata [1:], dtype=float)
xdata=data[:,0]
ydata=data[:,1]


with open ('params.txt', mode='r') as file:
	rawdata2=list(csv.reader(file, delimiter='\t'))
i=0
guess=[]
for i in range(4):
     k = np.float (rawdata2 [i] [1])
     guess.append(k)
     i=i+1

'''Use of the analytical equation that describes the TTOR model.'''
def fourPL2(x, I, b, Dc, a):
     return (I*(1-(np.real(lambertw(b* np.exp (b) * np.exp (-x/Dc)))/b)**a))


'''Fitting the experimental and the theoretical curves. The parameter bounds are set.'''
params, params_covariance = optimization.curve_fit(fourPL2, xdata, ydata, guess, bounds= (0, [np.inf, np.inf, np.inf, 0.999]), maxfev=5000)


#Output data
window=Tk()
window.eval('tk::PlaceWindow %s center' % window.winfo_toplevel())
window.withdraw()
messagebox.showinfo("Parameters", "The values of the fitting params are:\n I=%5.3f\n B=%5.3f\n Dc=%5.3f\n a=%5.3f. \n\n They have been saved in a .txt file in your folder." %tuple (params))

window.deiconify()
window.destroy()
window.quit()


'''FOM Calculation.'''
Ith=[]

for i in xdata:
	values=fourPL2(i, *params)
	Ith.append(values)
FOM = 100*(np.sum(abs(Ith-ydata))/(np.sum(ydata)))
print ('\n FOM=%0.3f' %FOM , '%')


'''The ascii file is written and includes the values of the theoretical curve'''
zip(xdata,Ith)
with open('Theoretical curve.txt', 'w', newline="") as f:
     writer = csv.writer(f, delimiter='\t')
     writer.writerow([rawdata [0] [0], "Ith (A.U.)"])
     writer.writerows(zip(xdata,Ith))

'''Creation of the txt file that includes the parameter values after the fitting.'''
parameters=['I=', 'B=', 'Dc=', 'a=']
with open('Fitting Params.txt', 'w', newline="") as f:
    writer = csv.writer(f, delimiter='\t')
    n=0
    for n in range(4):
         writer.writerow([parameters[n], params[n]])
         n=n+1    

'''Graphical plotting and data visualization.'''
x_min, x_max = np.amin(xdata), np.amax(xdata)
xs = np.linspace(x_min, x_max, 1000)
plt.scatter(xdata, ydata, label='experimental data \nFOM=%5.3f' %FOM)
plt.plot(xs, fourPL2(xs, *params),color='y', label='fit params:\n I=%5.3f\n B=%5.3f\n Dc=%5.3f\n a=%5.3f' %tuple (params))

plt.title ('DOSE RESPONSE')
plt.xlabel(np.str(rawdata [0] [0]))
plt.ylabel(rawdata [0] [1])
plt.legend()
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.savefig ('plot.png', bbox_inches='tight')
plt.show()


