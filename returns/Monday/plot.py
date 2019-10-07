import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import signal
file = 'LR04stack.csv'
time, val = [],[]
timed=[]
with open(file) as f:
    data = csv.reader(f,delimiter=',')
    for row in data:
        time.append(float(row[0]))
        val.append(float(row[1]))
#for times in time:
#    timed.append((times-5000)*1000)

plt.figure(2)
plt.grid()
power=(np.fft.fft(val))
overwhelming=np.abs(power)**2
plt.plot(overwhelming)

x=np.linspace(0,5320,5321)
vals=np.interp(x,time,val)
final=np.abs(np.fft.fft(vals))**2
y=np.fft.fftfreq(len(vals))
# THIS IS WHERE I GOT HELP
a,b=signal.periodogram(vals)
plt.figure(3)
plt.subplot(121)
plt.loglog(1/y,final)
plt.subplot(122)
plt.loglog(1/a,b)
#plt.scatter(x,vals,s=2)


plt.figure(1)
plt.scatter(x,vals,s=2)
plt.scatter(time, val, s=2)
plt.grid()
plt.ylim(5.2, 2.5)

plt.show()

