import numpy as np
import matplotlib.pyplot as plt
from math import exp, sin, pi, cos, log
import random
h=1.05e-27
c=3.e10
k=1.38e-16
Em,S=np.genfromtxt('S+sigma.txt',unpack=True) #выгружаем эффективную площадь
kev_to_erg = 1.6e-9
R=np.genfromtxt('radius.txt',unpack=True) # распределение по расстояниям
N=0 #число фотонов для данного объекта
Nph=np.zeros(620000) #массив значений энергии для logNlogS
q=np.zeros(620000)#массив для числа объектов
Nph[0]=10
f=open('N(S).txt','w')
for w in range(619999):
    Nph[w+1]=Nph[w]+10 #заполняем значения энергий
for i in range(10000): #основной цикл
    a=random.randint(0,8000) #определяем R
    R0=R[a]
    if R0<100*3.086e18: # определяем поглощение
        Nh=0.1*R0
    else:
        Nh=0.1*100*3.086e18+1*(R0-100*3.086e18)
    print('R0=',R0)
    b=random.random() #определяем массу, b - вероятность попадения в определененный диапазон
    #e=random.randint(10000,100000000) #определяем диапазон возрастов
    print('m=',b) 
    if b<0.04: #выбор файла с данными
        Age,r,Teff,L=np.genfromtxt('ResM04.txt',unpack=True)
    elif b<0.15:
        Age,r,Teff,L=np.genfromtxt('ResM05.txt',unpack=True)
    elif b<0.52:
        Age,r,Teff,L=np.genfromtxt('ResM06.txt',unpack=True)
    elif b<0.81:
        Age,r,Teff,L=np.genfromtxt('ResM07.txt',unpack=True) 
    elif b<0.91:
        Age,r,Teff,L=np.genfromtxt('ResM08.txt',unpack=True) 
    elif b<0.94:
        Age,r,Teff,L=np.genfromtxt('ResM09.txt',unpack=True) 
    else:
        Age,r,Teff,L=np.genfromtxt('ResM10.txt',unpack=True)
    e=random.randint(Age[0],Age[-1])
    u=0
    while Age[u]<1000000:
        u=u+1
    if e<1000000: # разыгрываем значение возраста из искомого диапазона после чего выбираем строку из файла для данного значения возраста
        d=random.randint(0,u-1)
    else:
        d=random.randint(u,len(Age)-1)
    print('teff=',Teff[d])
    r0=r[d]*10000 #переводим радиус в сантиметры
    for j in range(1012) :  #вылетает из-за слишком большой экспоненты  
            if Em[j]<0.284: #определяем коэффициенты для поглощения
                c0=34.6
                c1=267.9
                c2=-476.1
            elif Em[j]<0.4:
                c0=78.1
                c1=18.8
                c2=4.3
            elif Em[j]<0.532:
                c0=71.4
                c1=66.8
                c2=-51.4
            elif Em[j]<0.707:
                c0=95.5
                c1=145.8
                c2=-61.1
            elif Em[j]<0.867:
                c0=308.91
                c1=-380.6
                c2=294.
            elif Em[j]<1.303:
                c0=120.6
                c1=169.3
                c2=-47.7
            elif Em[j]<1.840:
                c0=141.3
                c1=146.8
                c2=-31.5
            elif Em[j]<2.471:
                c0=202.7
                c1=104.7
                c2=-17
            elif Em[j]<3.21:
                c0=342.7
                c1=18.7
                c2=0
            elif Em[j]<4.038:
                c0=352.2
                c1=18.7
                c2=0
            elif Em[j]<7.111:
                c0=433.9
                c1=-2.4
                c2=0.75
            elif Em[j]<8.331:
                c0=629
                c1=30.9
                c2=0
            else :
                c0=701.2
                c1=25.2
                c2=0
            sigma=(c2/Em[j]+c1/(Em[j]*Em[j])+c0/(Em[j]*Em[j]*Em[j]))/1e24*Nh    
            if ((Em[j]*kev_to_erg)/k/Teff[d])>680: # чтобы не вылетало из-за экспоненты
                break
            else:
                N=N+(2*r0**2*(Em[j]*kev_to_erg)**2/R0**2/h**3/c**2/(exp((Em[j]*kev_to_erg)/k/Teff[d])-1))*exp(-sigma)*S[j]*(Em[j+1]-Em[j])*kev_to_erg #считаем число фотонов
    print('N=',N*2500)
    for x in range(620000): #считаем logNlogS
        if N*2500>Nph[x]:
            q[x]=q[x]+1
        else:
            N=0
            break
    N=0
for z in range(620000): #вывод в файл
    f.write(str(Nph[z])+' ')
    f.write(str(q[z])+'\n')
f.close()
        
    

    
    
