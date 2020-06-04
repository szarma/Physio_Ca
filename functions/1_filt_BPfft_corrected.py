import numpy as np
from scipy.fftpack import rfft, irfft, fftfreq
import matplotlib.pyplot as plt
from sys import argv


lowBF=0.03
highBF=4.0


Tzac=1200  # zacetek vizualizacije v izhodnih slikah, v sekundah
Tkon=2400  # konec vizualizacije v izhodnih slikah, v sekundah

data=np.loadtxt(argv[1])
time   = data[:-1,0]  #za time uporabi celoten 0-ti stolpec minus zadnji zapis
signal = data[:-1,1:len(data[0])] #za signal pa uporabi vse podatke (-zadnji zapis), razen podatke v prvem stolpcu
# sampling=np.loadtxt("sampling.txt")
sampling = 1/np.diff(time).mean()


'''signal=data[:-1,:]
time=[]
for i in range(len(signal)):
    time.append(i/sampling)'''

Tzac1=int(Tzac*sampling)
Tkon1=int(Tkon*sampling)

#morebitni predhodnji linearni detrend
signalDT=np.zeros([len(signal),len(signal[0])])

'''for i in range(len(signal[0])): 
    koef=(signal[(len(signal))-1,i]-signal[1,i])/(len(signal)-2)
    for j in range(len(signal)):
        signalDT[j][i]=signal[j][i] - koef*(j)
datoteka=open("detrended.dat","w+")
for i in range(len(signalDT)):
    for j in range(len(signalDT[0])-1):
        signal[i][j]=signalDT[i][j]     # vnos detrendanega signala
        print>>datoteka,signalDT[i][j],    #   vejica -> presledek, no new line
    print>>datoteka,signalDT[i][len(signalDT[0])-1]     # brez vejice ->new line
datoteka.close()'''

data=[]
W=[] 
FS=[] 
for i in range(len(signal[0])):
    Wt = fftfreq(signal[:,i].size, d=time[1]-time[0])
    W.append(Wt) #i-temu zapisu v listi dodelima vrednosti Wt
    f_signal = rfft(signal[:,i])
    
    cut_f_signal = f_signal.copy()
    cut_f_signal[(Wt<lowBF)] = 0 
    cut_f_signal[(Wt>highBF)] = 0
    cut_signal = irfft(cut_f_signal)
    FS.append(cut_signal)
    print ('Obdelujem TS: ',i)
       
    plt.subplot(211)
    plt.plot(time,signal[:,i])
    plt.subplot(212)
    plt.plot(time,cut_signal)
    plt.show()
    plt.close()
    
    plt.subplot(311)
    plt.plot(time,signal[:,i],linewidth=0.5,color='dimgrey')
    plt.title(r'trace $ %ld $'%(float(i)/1.0),fontsize=8)
    plt.rc('font', size=6)
    plt.subplot(312)
    miny=np.min(signal[(Tzac1):(Tkon1),i])
    maxy=np.max(signal[(Tzac1):(Tkon1),i])
    plt.plot(time,signal[:,i],linewidth=0.67,color='dimgrey')
    plt.xlim([Tzac,Tkon])
    plt.ylim([miny-0.05,maxy+0.05])
    plt.subplot(313)
    miny=np.min(cut_signal[Tzac1:Tkon1])
    maxy=np.max(cut_signal[Tzac1:Tkon1])
    plt.plot(time,cut_signal,linewidth=0.67,color='dimgrey')
    plt.xlim([Tzac,Tkon])
    plt.ylim([miny*1.05-0.05,maxy*1.05+0.05])
    plt.xlabel('time (s)', fontsize=6)
    plt.savefig("filt_traces\\filt_trace_%d.jpg"%(i+0),dpi=200,bbox_inches = 'tight')
    plt.close()  

zanemari=[]

datoteka=open("DATAfilt.txt","w+")
for i in range(len(signal)):
    for j in range(len(signal[0])-1):
        if j not in zanemari: ## dodano za zanemari
            print(datoteka,FS[j][i],end=" " )    #   vejica -> presledek, no new line
    print(datoteka,FS[len(signal[0])-1][i])     # brez vejice ->new line
datoteka.close()
