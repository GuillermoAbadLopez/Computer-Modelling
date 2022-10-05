import numpy as np 
import matplotlib.pylab as plt
import matplotlib.animation as animation
from decimal import Decimal, ROUND_DOWN

#FOTOS: per comprobar matrius
def heatmap2d(arr: np.ndarray):
    plt.imshow(arr, cmap='magma',vmin=0, vmax=1, animated=True)
    plt.colorbar()
    plt.show()

#DADES BÀSIQUES
L = 20 #20x20m 
#ample carretera = 3m
dx, dy, dt = 0.1 , 0.1, 0.09 #0.1m; 0.09s 
X, Y = np.int(L/dx+1) , np.int(L/dy+1)
Tmax = 1000
e = 3/dx #gruix ha de ser impar, ample carretera
ee = (3-0.5)/dx #gruix diferent per a que la rotonda no divergeixi


####################################################################################################################################################################################################################################
#INTERSECCIO:
####################################################################################################################################################################################################################################

#ZONAS PERMITIDAS
ZOi = np.zeros((X,X), int)
ANCHOFINALi=1/dx       #Zona de carril que actua como sumideros de U
ANCHOINICIALi=5/dx     #Zona de carril que actua como fuente de U

for i in range(X):
    for j in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if i>ANCHOFINALi-1: #sumidero en salidas, evitamos colapso
            ZOi[i,j]=1
for j in range(X):
    for i in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if j<(X-ANCHOFINALi): #sumidero en salidas, evitamos colapso
            ZOi[i,j]=1
#heatmap2d(ZOi)

#PESO DIRECCIONAL PARA V EN Y
VYi = np.zeros((X,X))
for i in range(X):
    for j in range(np.int(X/2-e/2),np.int(X/2+e/2)): #ancho carril
        if i<np.int(X/2-e/2) or i > np.int(X/2+e/2):
            VYi[i,j] = -1 #puesto que tiene que ir hacia arriba, i nuestro eje esta definido positivo hacia abajo
        elif j<np.int(X/2+e/2+1) and j>np.int(X/2+e/4) and i<np.int((X-1)-j) and i<np.int(X/2):                           #Per evitar les colisions amb la esquina inlino la matriu de velocitats
            VYi[i,j] = -1 #idem "
        elif j<np.int(X/2+e/2+1) and j>np.int(X/2+e/4) and i>np.int((X-1)-j) and i<np.int(X/2-e/4):
            VYi[i,j] = 0 #para evitar divergencia en el cruce
        else:
            VYi[i,j] = -1/np.sqrt(2)-0.1 #peso en el cruce
#heatmap2d(VYi)
            
#PESO DIRECCIONAL PARA V EN X (idem ")
VXi= np.zeros((X,X))
for j in range(X):
    for i in range(np.int(X/2-e/2),np.int(X/2+e/2)):
        if j<np.int(X/2-e/2) or j > np.int(X/2+e/2):
            VXi[i,j] = 1
        elif j<np.int(X/2+e/2+1) and j>np.int(X/2+e/4) and i<np.int((X-1)-j) and i<np.int(X/2):
            VXi[i,j] = 0                                                                                                    #Per evitar les colisions amb la esquina inlino la matriu de velocitats
        elif j<np.int(X/2+e/2+1) and j>np.int(X/2+e/4) and i>np.int((X-1)-j) and i<np.int(X/2-e/4):
            VXi[i,j] = 1
        else:
            VXi[i,j] = 1/np.sqrt(2)-0.1
#heatmap2d(VXi)
            
#CCII
CCIItotespaii=0.0      #Condicions inicials a tot el espai
CCIIfontsix=0.4        #Condicions inicials de la font carril horitzonal
CCIIfontsiy=0.4        #Condicions inicials de la font carril vertical
CCFFfontsi=0.4         #Condicions frontera de les fonts

#METODE
dg = dt/dx
Ui = np.zeros((X,X)) #ccii
Uinew = np.zeros((X,X)) #iterativa
UTi= np.zeros((Tmax,X,X)) #tensor

#Aplicamos ccii
for i in range(X):#modificadas para simul iba X
    for j in range(Y):
        Ui[i,j] = ZOi[i,j]*CCIItotespaii #coxtes a tot l'espai a t=0        
        if j<ANCHOINICIALi: #cotxes a les fonts a t=0
            Ui[i,j] = ZOi[i,j]*CCIIfontsix
        if i>X-2-ANCHOINICIALi:
            Ui[i,j] = ZOi[i,j]*CCIIfontsiy
#heatmap2d(Ui)
            
for t in range(Tmax):
    for i in range(1,X-1):#no toquem bordes 
        for j in range(1,Y-1): #no toquem bordes, es queden tots com Uinew incial, es a dir 0's
            Uinew[i,j] = Ui[i,j]-dg/4*VYi[i,j]*(Ui[i+1,j]**2-Ui[i-1,j]**2)+((VYi[i,j]*dg)**2)/8*((Ui[i,j]+Ui[i+1,j])*(Ui[i+1,j]**2-Ui[i,j]**2)-(Ui[i,j]+Ui[i-1,j])*(Ui[i,j]**2-Ui[i-1,j]**2))-dg/4*VXi[i,j]*(Ui[i,j+1]**2-Ui[i,j-1]**2)+((VXi[i,j]*dg)**2)/8*((Ui[i,j]+Ui[i,j+1])*(Ui[i,j+1]**2-Ui[i,j]**2)-(Ui[i,j]+Ui[i,j-1])*(Ui[i,j]**2-Ui[i,j-1]**2))

    for i in range(X):
        for j in range(Y):
            Ui[i,j] = Uinew[i,j]
            """                              #Posar cometes aquí, si no es vol terme de fuente (paquete)
            if j<ANCHOINICIALi: 
                Ui[i,j] = ZOi[i,j]*CCFFfontsi
            if i>X-2-ANCHOINICIALi:
                Ui[i,j] = ZOi[i,j]*CCFFfontsi 
            """                              #Acabar les cometes aquí
            
            
            
            
    #PRINTEO RÁPID
    IntervaloPrint=1
    to=t%IntervaloPrint
    if to==0:
        heatmap2d(Ui)



"""
UTi[t]=Ui
#VELOCITATS
VTi= np.zeros((Tmax,X,X)) #tensor
for t in range(Tmax):
            VTi[t]=(1-(((1-UTi[t])*(135/2))/135))*30-15
"""            
            
"""
#ANIMACIÓ   
fig = plt.figure()
ims = []
for i in range(Tmax):
    im = plt.imshow(VTi[i], cmap='magma',vmin=0, animated=True)
    ims.append([im])
    
plt.colorbar()
ani = animation.ArtistAnimation(fig, ims, blit=True)

ani.save('INTER.mp4')

plt.show()
"""        





   
####################################################################################################################################################################################################################################
#ROTONDA :
####################################################################################################################################################################################################################################

r = (2+3)/dx #radi par (rotonda+carretera)
a = X/2-0.5 #centrada a (a,a)
x, y = np.indices((X, X))
ZOr = np.array(np.abs(np.hypot(a - x, a - y)-r) < e/2).astype(int)
ZOrr = np.array(np.abs(np.hypot(a - x, a - y)-r) < (e)/2).astype(int)
ZOrrr = np.array(np.abs(np.hypot(a - x, a - y)-r) < (ee)/2).astype(int)
ANCHOFINALr=1/dx         #Llargaria dels sumideros on desapareixen les rhos
ANCHOINICIALr=3/dx       #Llargaria de les fonts on es creen les rhos

#zones permeses y CCII
for i in range(X):
    for j in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if i< np.int(X/2-0.5-r-2) and i>(ANCHOFINALr-1):
            ZOr[i,j]=1
        if i> np.int(X/2+0.5+r+1):
            ZOr[i,j]=1
for j in range(X):
    for i in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if j< np.int(X/2-0.5-r-2):
            ZOr[i,j]=1
        if j > np.int(X/2+0.5+r+1) and j<(X-ANCHOFINALr):
            ZOr[i,j]=1
            
#zones permeses2.0
for i in range(X):
    for j in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if i< np.int(X/2-0.5-r-2) and i>(ANCHOFINALr-1):
            ZOrrr[i,j]=1
        if i> np.int(X/2+0.5+r+1):
            ZOrrr[i,j]=1
for j in range(X):
    for i in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if j< np.int(X/2-0.5-r-2):
            ZOrrr[i,j]=1
        if j > np.int(X/2+0.5+r+1) and j<(X-ANCHOFINALr):
            ZOrrr[i,j]=1

                    
#MATRIU D'ANGLES: i li podem aplicar dp cos() i sin() com a funcions a cada angle i per cada matriu particular q volguem
s=X #=X
p = np.zeros((s,s))
q = np.zeros((s,s))
RADIAL=0.1   #El flow s'en va cap a dins de la rotonda en augmentar aquest número

for i in range(s): #()
    for j in range(s):
        q[j,i]= -(i-(s/2-0.499999))/(j-(s/2-0.499999))
        p[i,j] = np.arctan2(1,1/q[j,i])
        if i>s/2-0.5:
            p[i,j] = p[i,j] + np.pi
   
     
#VX 
VXr=np.zeros((s,s))
Carrx=1.5/dx
PadentroElCarril=2
#Cos/Sin(Angles)* rotonda
for i in range(s): #()
    for j in range(s):
      """   """                                         #CAMBIAR COMETES AMUNT O ABAIX DEPENEN DE SI VOLS LA MEITAT O TOT
      if j<i:
          VXr[i,j]=ZOrr[i,j]*np.sin(p[i,j]+RADIAL)*(-1)
      else:
          VXr[i,j]=ZOrr[i,j]*np.sin(p[i,j])*(-1)
      """     """                                        #CAMBIAR COMETES AMUNT O ABAIX DEPENEN DE SI VOLS LA MEITAT O TOT
      VXr[i,j]=ZOrr[i,j]*np.sin(p[i,j]+RADIAL)*(-1)
      
            
#(Cos/Sin(Angles)* rotonda)+carriles                                    MIRAR LIMITES CARRILES!!!!!
for j in range(X):
    for i in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if j< np.int(X/2-0.5-r-e/2+Carrx*PadentroElCarril) and j<np.int(X/2-5-i/2):   #inclino els carrils a la sortida
            VXr[i,j]=1
        if j > np.int(X/2+0.5+r+e/2-Carrx-1)and j>np.int(0.6*X+i*0.36):
            VXr[i,j]=1
for i in range(X):
    for j in range(np.int(X/2-e/2),np.int(X/2+e/2)):
        if i< np.int(X/2-0.5-r-e/2+Carrx*PadentroElCarril) and i<np.int(X/2-5-j/2) and i<np.int(0*X-5+j/2):
            VXr[i,j]=0
        if i> np.int(X/2+0.5+r+e/2-Carrx-1) and i>np.int(0.6*X+j*0.35):
            VXr[i,j]=0

                
#VY
VYr=np.zeros((s,s))
Carry=1.5/dy
#Cos/Sin(Angles)* rotonda
for i in range(s): #()
    for j in range(s):
        """    """                     #CAMBIAR COMETES AMUNT O ABAIX DEPENEN DE SI VOLS LA MEITAT O TOT
        if j<i:
            VYr[i,j]=-ZOrr[i,j]*np.cos(p[i,j]+RADIAL)
        else:
            VYr[i,j]=-ZOrr[i,j]*np.cos(p[i,j])     
        """    """                      #CAMBIAR COMETES AMUNT O ABAIX DEPENEN DE SI VOLS LA MEITAT O TOT
        VYr[i,j]=-ZOrr[i,j]*np.cos(p[i,j]+RADIAL)
        
#(Cos/Sin(Angles)* rotonda)+carriles                                    MIRAR LIMITES CARRILES!!!!!
for i in range(X):
    for j in range(np.int(X/2-e/2),np.int(X/2+e/2)): #amplada carril 
        if i< np.int(X/2-0.5-r-e/2+Carry*PadentroElCarril) and i<np.int(X/2-5-j/2) and i<np.int(0*X-5+j/2):
            VYr[i,j]=-1
        if i> np.int(X/2+0.5+r+e/2-Carry-1) and i>np.int(0.6*X+j*0.35):
            VYr[i,j]=-1
for j in range(X):
    for i in range(np.int(X/2-e/2),np.int(X/2+e/2)):
        if j< np.int(X/2-0.5-r-e/2+Carry*PadentroElCarril) and j<np.int(X/2-5-i/2):   #inclino els carrils a la sortida
            VYr[i,j]=0
        if j > np.int(X/2+0.5+r+e/2-Carry-1) and j>np.int(0.6*X+i*0.36):
            VYr[i,j]=0

            

#CCII
CCIItotespair=0.0      #Condicions inicials a tot el espai
CCIIfontsr=0.7        #Condicions inicials de les fonts
CCFFfontsr=0.7        #Condicions frontera de les fonts


#METODE:
dg = dt/dx
Tmax = 1500
#CCII
Ur=np.zeros((X,X))
Urnew=np.zeros((X,X))
UTr= np.zeros((Tmax,X,X)) #tensor
for i in range(X):
    for j in range(Y):
        Ur[i,j] = ZOr[i,j]*CCIItotespair
        if j<ANCHOINICIALr: 
            Ur[i,j] = ZOr[i,j]*CCIIfontsr
        if i>X-2-ANCHOINICIALr:
            Ur[i,j] = ZOr[i,j]*CCIIfontsr
            
for t in range(Tmax):
    for i in range(1,X-1):
        for j in range(1,Y-1):
            Urnew[i,j] = Ur[i,j]-dg/4*VYr[i,j]*(Ur[i+1,j]**2-Ur[i-1,j]**2)+((VYr[i,j]*dg)**2)/8*((Ur[i,j]+Ur[i+1,j])*(Ur[i+1,j]**2-Ur[i,j]**2)-(Ur[i,j]+Ur[i-1,j])*(Ur[i,j]**2-Ur[i-1,j]**2))-dg/4*VXr[i,j]*(Ur[i,j+1]**2-Ur[i,j-1]**2)+((VXr[i,j]*dg)**2)/8*((Ur[i,j]+Ur[i,j+1])*(Ur[i,j+1]**2-Ur[i,j]**2)-(Ur[i,j]+Ur[i,j-1])*(Ur[i,j]**2-Ur[i,j-1]**2))
    for i in range(X):
        for j in range(Y):
            Ur[i,j] = Urnew[i,j]*ZOr[i,j]
    for i in range(X):
        for j in range(Y):
            if Ur[i,j] < 0:
                Ur[i,j]=-Ur[i,j]
            elif Ur[i,j] > 0.95:
                Ur[i,j]=0.9

                                        #Possar cometes aquí si no volem fuente
            if j<ANCHOINICIALr: 
                Ur[i,j] = ZOr[i,j]*CCFFfontsr
            if i>X-2-ANCHOINICIALr:
                Ur[i,j] = ZOr[i,j]*CCFFfontsr
                                        #Possar cometes aquí si no volem fuente
    
    
    
    #PRINTEO RÁPIDO
    IntervaloPrint=1
    to=t%IntervaloPrint
    if to==0:
        heatmap2d(Ur)
        
        
    UTr[t]=Ur    
"""      
#ANIMACIÓ   
fig = plt.figure()
ims = []
for i in range(Tmax):
    im = plt.imshow(UTi[i], cmap='magma',vmin=0, vmax=1, animated=True)
    ims.append([im])
    
plt.colorbar()
ani = animation.ArtistAnimation(fig, ims, interval=10, blit=True)

ani.save('INTER.mp4')

plt.show()
"""