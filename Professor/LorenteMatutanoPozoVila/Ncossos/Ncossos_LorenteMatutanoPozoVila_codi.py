# -*- coding: utf-8 -*-
"""MNII-Sistema Caótico

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/14iiP30BWLiEF26TD1mhRqgXp_OnaYC4R

Preàmbul
"""

import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size': 12})
import matplotlib.pyplot as plt
from google.colab import files

from google.colab import drive
drive.mount('/content/gdrive')

def unit(x,y):                                      #Vector unitari dels angles definits
  return(np.array([np.sin(x)*np.cos(y),np.sin(x)*np.sin(y),-np.cos(x)]))
def mod(x):                                         #Módul d'un vector
  return(np.sqrt(sum(x*x)))
def f_1(γ,α,z,r,r_1,r_2):                           #Acceleració d'm_1
  return(γ*(1-α/r)*(r_2-np.dot(r_1,r_2)*r_1)+np.dot(r_1,z)*r_1-z-(mod(v_1)**2)*r_1)
def f_2(χ,γ,α,z,r,r_1,r_2):                         #Acceleració d'm_2
  return(-(χ*γ*(1-α/r)*(r_2-r_1)+z))
z=np.array([0,0,1])     #Vector unitari z, simplificarà els càlculs

def rep(R_1,R_2,q,size,t):                          #Funció representació gràfica en 3D
  R_1=np.transpose(R_1)
  R_2=np.transpose(R_2)

  fig = plt.figure()
  ax = plt.axes(projection='3d')
  ax.scatter3D(R_1[0], R_1[1], R_1[2],color='indigo',s=5*np.linspace(0,1,len(R_1[0])))  #Posició i estela d'm_1
  ax.scatter3D(R_2[0], R_2[1], R_2[2],color='violet',s=5*np.linspace(0,1,len(R_2[0])))  #Posició i estela d'm_2
  ax.scatter3D([0],[0],[0],color='white',edgecolor='black',s=20,marker='v')             #Origen
  ax.plot3D([0,q[0]],[0,q[1]],[0,q[2]],color='black')                             #Barra rígida
  ax.plot3D([q[0],q[3]],[q[1],q[4]],[q[2],q[5]],color='black', dashes=[2*r,2*r])#Molla
  #Format del gràfic:
  ax.set_xlim3d(-size,size)
  ax.set_ylim3d(-size,size)
  ax.set_zlim3d(-size,size)
  ax.set_xlabel('$x/L$')
  ax.set_ylabel('$y/L$')
  ax.set_zlabel('$z/L$')
  ax.text(size, -size, 5*size/6, '$t=$'+str(round(t,1))+r'$\sqrt{L/g}$', color='black')
  ax.text(size, -size, size/2, r'$|r_1|=$'+str(round(mod(q[0:3]),2))+r'$L$', color='black')
  ax.view_init(30,45)
    #Generem i guardem la imatge:
  plt.savefig(str(round(t,1))+'.png')
  plt.savefig(f"{'/content/gdrive/MyDrive/PenduloMuelleTemp'}/"+'frame'+str(int(10*round(t,1)))+".png")

  R_1=np.transpose(R_1)
  R_2=np.transpose(R_2)

z=np.array([0,0,1])     #Definim el vector unitari z, simplificarà els càlculs.

"""Simulació"""

#Tots els ANGLES s'expressaran en radians, les LONGITUDS en L, TEMPS en sqrt(L/g) i VELOCITATS en sqrt(Lg)

#Propietats del sistema:

α=0.5                   #Valor l/L, que determina la mida del pèndol
γ=5                     #Valor k*L/g*m_1, que determina la força relativa de la molla
χ=10                    #Valor m_1/m_2, que determina quina de les masses es mourà més erràticament
dt=0.01                 #Pas de temps
tmax=10                 #Temps màxim que volem estudiar

#Condicions inicials:

#Les posicions hem considerat que és més senzill introduir-les en COORDENADES GENERALITZADES
θ_1=0
θ_2=0.2
φ_1=0.2
φ_2=0.2
r=1.2*α                 #Longitud inicial de la molla (no la natural)

#Les velocitats, però directament en COORDENADES CARTESIANES
v_1=np.array([0,0,0])
v_2=np.array([0,0,0])
w=np.append(v_1,v_2)    #Vector unificat de velocitats

t=0.

#Corresponents vectors cartesians:

r_1=unit(φ_1,θ_1)
r_2=unit(φ_1,θ_1)+r*unit(φ_2,θ_2)
q=np.append(r_1,r_2)    #Vector unificat de posicions

#Propietats de la representació:

GIF=True                #Volem fer GIF? Sí=>TRUE // No=>FALSE
N=50                    #Número de punts passats que volem representar en cada fotograma.
f=0.1                   #Període que volem entre fotogrames
F=int(f/dt)             #Ho expressem en numero de passes temporals entre fotogrames
count=0                 #Conteig de passes temporals

R_1=np.array([q[0:3]])  #En aquestes matrius guardarem les posicions a representar
R_2=np.array([q[3:6]])

if GIF==True:           #Representem l'instant inicial
  rep(R_1,R_2,q,3,t)

e_0=mod(w[0:3])**2/2+q[2]+(mod(w[3:6])**2/2+q[5])/χ+(γ/2)*(r-α)**2
T=np.array(0)           #Matrius de temps i energies
E=np.array(e_0)

while t<=tmax:
  #Fem el càlcul numèric
  dq_1=w*dt
  dw_1=np.append(f_1(γ,α,z,mod(q[3:6]-q[0:3]),q[0:3],q[3:6]),f_2(χ,γ,α,z,mod(q[3:6]-q[0:3]),q[0:3],q[3:6]))*dt
  
  q_aux=q+dq_1/2
  dq_2=(w+dw_1/2)*dt
  dw_2=np.append(f_1(γ,α,z,mod(q_aux[3:6]-q_aux[0:3]),q_aux[0:3],q_aux[3:6]),f_2(χ,γ,α,z,mod(q_aux[3:6]-q_aux[0:3]),q_aux[0:3],q_aux[3:6]))*dt
  
  q_aux=q+dq_2/2
  dq_3=(w+dw_2/2)*dt
  dw_3=np.append(f_1(γ,α,z,mod(q_aux[3:6]-q_aux[0:3]),q_aux[0:3],q_aux[3:6]),f_2(χ,γ,α,z,mod(q_aux[3:6]-q_aux[0:3]),q_aux[0:3],q_aux[3:6]))*dt
  
  q_aux=q+dq_3
  dq_4=(w*dw_3)*dt
  dw_4=np.append(f_1(γ,α,z,mod(q_aux[3:6]-q_aux[0:3]),q_aux[0:3],q_aux[3:6]),f_2(χ,γ,α,z,mod(q_aux[3:6]-q_aux[0:3]),q_aux[0:3],q_aux[3:6]))*dt

  q=q+(dq_1+2*dq_2+2*dq_3+dq_4)/6
  w=w+(dw_1+2*dw_2+2*dw_3+dw_4)/6

  t=t+dt
  count+=1              #Anotem que hem fet un pas de temps

  T=np.append(T,t)
  E5=np.append(E5,mod(w[0:3])**2/2+q[2]+(mod(w[3:6])**2/2+q[5])/χ+(γ/2)*(mod(q[3:6]-q[0:3])-α)**2)

  if count==F:          #Si ha passat el temps dessitjat entre fotogrames:
    if GIF==True:
      #Anotem la nova informació
      R_1=np.append(R_1,[q[0:3]],axis=0)
      R_2=np.append(R_2,[q[3:6]],axis=0)
      #I la representem
      rep(R_1,R_2,q,3,t)
      #Eliminem els punts de les matrius que ja no farem servir
      if len(R_1)==N:
        R_1=R_1[1:]
        R_2=R_2[1:]

    count=0             #Restablim el conteig de passes

"""Gràfic energies"""

plt.plot(T,E5/e_05,color='darkred'   ,label='φ=0.5')
plt.plot(T,E4/e_04,color='coral'     ,label='φ=0.4')
plt.plot(T,E3/e_03,color='darkorange',label='φ=0.3')
plt.plot(T,E2/e_02,color='gold'      ,label='φ=0.2')
plt.plot(T,E1/e_01,color='limegreen' ,label='φ=0.1')
plt.plot(T,E0/e_00,color='darkgreen' ,label='φ=0.0')
plt.xlabel(r'$t$ $\left(\sqrt{L/g}\right)$')
plt.ylabel(r'$E(t)$ / $E(0)$')
plt.xlim(0,50)
plt.legend()
plt.show()