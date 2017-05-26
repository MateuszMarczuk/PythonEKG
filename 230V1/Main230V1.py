import wfdb
import csv
import pywt
import matplotlib.pyplot as pyplot
import sys
import numpy as np
import sys


hea = ""
if len(sys.argv) > 1:
    hea = str(sys.argv[1])

record = wfdb.rdsamp('230')
sig, fields = wfdb.srdsamp('230', channels = [0])
annotation = wfdb.rdann('230', 'atr')
# signalMLII = record.p_signals[:, 0]
#         signalV5 = record.p_signals[:, 1]
#         rate = record.fs;
#         t = np.arange(0, len(signalMLII) / rate, 1 / rate)
#         qrs_ref = ann.annsamp[1:] / rate
wfdb.plotrec(record, annotation = annotation, title='Record 100 from MIT-BIH Arrhythmia Database', timeunits = 'seconds')

Tree = pywt.wavedec(data=sig, wavelet='db2', level=2, mode='symmetric', axis=-2)

for i in range(0, len(Tree[0])):
   if (Tree[0][i] <0.4):
      Tree[0][i]=0
newTree = [Tree[0], Tree[1]*0, Tree[2]*0]


TreeHaar = pywt.wavedec(data=sig, wavelet='dmey', level=2, mode='symmetric', axis=-2)

for i in range(0, len(TreeHaar[0])):
   if TreeHaar[0][i] < 0.4:
       TreeHaar[0][i]=0



newTreedmey = [TreeHaar[0], TreeHaar[1]*0, TreeHaar[2]*0]


Tree3 = pywt.wavedec(data=sig, wavelet='coif1', level=3, mode='symmetric', axis=-2)

for i in range(0, len(Tree3[0])):
   if Tree3[0][i] < 0.4:
      Tree3[0][i]=0

newTreeCoif = [Tree3[0], Tree3[1]*0, Tree3[2]*0, Tree3[3]*0]




recSignal=pywt.waverec(newTree,'db2',axis=-2)
recSignalHaar=pywt.waverec(newTreedmey,'dmey',axis=-2)
recSignal3=pywt.waverec(newTreeCoif,'coif1',axis=-2)


pyplot.figure(1)
pyplot.plot(recSignal[:2000])
pyplot.title("Falka db2, level = 5")
pyplot.show()
pyplot.figure(2)
pyplot.plot(recSignalHaar[:2000])
pyplot.title("Falka dmey, level = 2")
pyplot.show()
pyplot.figure(3)
pyplot.plot(recSignal3[:2000])
pyplot.title("Falka coif2, level = 5")
pyplot.show()

def maximum (recSignal):
    max=[]
    for i in range(0, len(sig) - 1):
        if (recSignal[i] > recSignal[i - 1]) and (recSignal[i] > recSignal[i + 1]):
            if (len(max) == 0) or (i - max[-1] > 0.4 * record.fs):
                max.append(i)
            elif recSignal[i] > recSignal[max[-1]]:
                max[-1] = i
    return max

def bp (max):
    bpm = []
    result = []
    for i in range(0, len(max)-1):
        liczba=max[i+1]-max[i]
        jednostka=liczba/360
        wynik=60/jednostka
        bpm.append(wynik)
        result.append([str(round(max[i])), str(round(bpm[i], 2))])

        print("Lokalizacja tonu =", round(max[i]), "|", "BMP =", round(bpm[i], 2))
    return result


def adno (max, probka=150):
    TP=0
    FP=0
    #FN =0
    for i in range(0 ,len(annotation.annsamp)): # po wszsytkich adnotacjach wejściowych
        adnotacja_wej=annotation.annsamp[i];
        pierwsza_trafiona=0;
        for j in range(0, len(max)):
            adnotacja_moja=max[j];
            if (adnotacja_wej - probka < adnotacja_moja) and (adnotacja_wej + probka > adnotacja_moja):
                if (pierwsza_trafiona==0): # jeżeli adnotacja jest pierwsza kolejna już nie będzie
                    pierwsza_trafiona = 1;
                    TP+=1;

                else: # jeżeli adnotacja nie jest pierwsza
                    FP+=1;

    return TP, FP

def FalseNegative(max):
    TP, FP = adno(max)
    FN=len(annotation.annsamp)- TP
    return FN

def Se (max):
    TP, FP = adno(max)
    SE = TP / (TP + FalseNegative(max))
    return SE

def RR (max):
    TP, FP = adno(max)
    PPV = TP / (TP + FP)
    return RR



def Suma (max, probka=150):
    SUMA=0
    TP, FP = adno(max)
    for i in range(0 ,len(annotation.annsamp)):
        adnotacja_wej=annotation.annsamp[i];
        for j in range(0, len(max)):
            adnotacja_moja=max[j];
            if (adnotacja_wej - probka < adnotacja_moja) and (adnotacja_wej + probka > adnotacja_moja):
                SUMA+=abs(adnotacja_moja-adnotacja_wej)/TP
    return SUMA




plikcsv="DB2.MLII.BPM.230.csv"
with open(plikcsv, 'w', newline='') as file:
    x = csv.writer(file, delimiter=';')
    x.writerows(bp(maximum(recSignal)))
    file.close()


plikcsv="DMEY.MLII.BPM.230.csv"
with open(plikcsv, 'w', newline='') as file:
    x = csv.writer(file, delimiter=';')
    x.writerows(bp(maximum(recSignalHaar)))
    file.close()


plikcsv="COIF1.MLII.BPM.230.csv"
with open(plikcsv, 'w', newline='') as file:
    x = csv.writer(file, delimiter=';')
    x.writerows(bp(maximum(recSignal3)))
    file.close()

def wynik (max):
    TP, FP = adno(max)
    FN=FalseNegative(max)
    SE = TP / (TP + FalseNegative(max))
    PPV = TP / (TP + FP)
    SUMA=Suma(max)
    csvwynik2=[]
    csvwynik2.append([str(TP),str(FP),str(FN),str(SE),str(PPV),str(SUMA)])
    return csvwynik2

print(wynik(maximum(recSignal)))
print(wynik(maximum(recSignalHaar)))
print(wynik(maximum(recSignal3)))


plikcsv2="Wynik230.DB2.MLII.csv"
with open(plikcsv2, 'w', newline='') as file:
    x = csv.writer(file, delimiter=';')
    x.writerows(wynik(maximum(recSignal)))
    file.close()

plikcsv2="Wynik230.DMEY.MLII.csv"
with open(plikcsv2, 'w', newline='') as file:
    x = csv.writer(file, delimiter=';')
    x.writerows(wynik(maximum(recSignalHaar)))
    file.close()

plikcsv2="Wynik230.COIF1.MLII.csv"
with open(plikcsv2, 'w', newline='') as file:
    x = csv.writer(file, delimiter=';')
    x.writerows(wynik(maximum(recSignal3)))
    file.close()



