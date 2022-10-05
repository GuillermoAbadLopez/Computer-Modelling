import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Polygon, Circle

#CONSTANTS
m = [1,100000000]             #llista de masses normalitzades
R = [5,10000]                       #llista de radis de les masses
c = 1.2                         #constant d'esmorteïment
w = 1                           #velocitat de l'aigua (normalitzada)

V = c/w                         #constant d'esmorteïment (normalitzada)
q0 = 1                          #alçada plana del mar (normalitzada)

k_q = 500000                     #coeficient de trasbalsament (moviment de marees d'un punt al següent)


#DISCRETITZACIÓ
T = 7000                        #nombre de Δt's
Z = 100                         #nombre de Δθ's
dt = 0.1                        #Δt
dz = 2*np.pi/Z                  #Δθ
gamma = 1.1*dt                  #γ = Δt/Δθ
D = 2                           #nombre de dimensions

    #N cossos
r = np.zeros((len(m),T,D))      #posició [massa][temps][eix]
v = np.zeros((len(m),T,D))      #velocitat [massa][temps][eix]

    #Marees
q = np.zeros((Z,T))             #alçada del mar [angle][temps]

t = dt*T/w                      #temps final
print('\n La simulació es duu a terme entre t=0 i t=%.2f (en unitats de temps normalitzades).\n' %t)


#CONDICIONS INICIALS
for i in range(Z):                                      #Posicions inicials del mar
    q[i][0] = q0
    q[i][1] = q0
                                                        #Posicions inicials/velocitats inicials
                                                        #dels N cossos
    #Massa 0
r[0][0][0] = 100000
r[0][0][1] = 0
v[0][0][0] = 0
v[0][0][1] = 50


    #Massa 2
r[1][0][0] = 0
r[1][0][1] = 0
v[1][0][0] = 0
v[1][0][1] = 0






#SIMULACIÓ
def RK4(k):                         #0, 1 o 1/2, depenent del vector del mètode RK4
    if (k==0) or (k==3):            #que multiplicarà els coeficients de la posició
        return k/3                  #i velocitat actualitzades
    else:
        return 1/2


def G(i,j,d):                                                           #atracció gravitatòria que experimenta
    f = 0                                                               #el cos i, en l'instant j, en l'eix d
    for l in [n for n in range(len(m)) if n != i]:
        div = 0
        for s in range(D):
            div += pow(r[i][j][s] - r[l][j][s],2)
        f += - m[l] * (r[i][j][d] - r[l][j][d]) / pow(div,1.5)
    return f


def g(i,j,k,B,d):                                                   #atracció gravitatòria adaptada al
    f = 0                                                           #mètode RK4 (inclou els increments
    dr = np.zeros((len(m),D))                                       #de posició en la matriu B)
    for l in range(len(m)):
        for s in range(d):
            dr[l][s] = dt*RK4(k)*B[l][s]

    for l in [n for n in range(len(m)) if n != i]:
        div = 0
        for s in range(D):
            div += pow(r[i][j][s] + dr[i][s] - r[l][j][s] - dr[l][s],2)
        f += - m[l] * (r[i][j][d] + dr[i][d] - r[l][j][d] - dr[l][d]) / pow(div,1.5)
    return f


for j in range(0,T-1):                          #càlcul de la posició del cos i en l'instant j+1
    C = np.zeros((4,len(m),2*D))
    for i in range(len(m)):                     #coeficient k del mètode RK4 per al cos i en cada eix
        for d in range(D):                      #C_k[massa][eix*2]
            C[0][i][d] = v[i][j][d]
            C[0][i][d+D] = G(i,j,d)             #primer coeficient

    for i in range(len(m)):
        for d in range(D):
            for k in range(1,4,1):
                C[k][i][d] = v[i][j][d] + RK4(k)*dt*C[k-1][i][d]    #següents coeficients
                C[k][i][d+D] = g(i,j,k,C[k-1,:,:],d)

        for d in range(D):                                          #actualització de posicions i velocitats
            v[i][j+1][d] = v[i][j][d] + (1/6)*dt*(C[0][i][d+D]+2*C[1][i][d+D]+2*C[2][i][d+D]+C[3][i][d+D])
            r[i][j+1][d] = r[i][j][d] + (1/6)*dt*(C[0][i][d]+2*C[1][i][d]+2*C[2][i][d]+C[3][i][d])




def Phi(l,j):               #força que experimenta el punt del mar l en l'instant j
    r_mar = [r[0][j][0] +(R[0]+q[l][j])* np.cos(2 * np.pi * l / Z), r[0][j][1] + (R[0]+q[l][j])* np.sin(2 * np.pi * l / Z)]
    U = [0,0]
    U_mar = [0,0]
    for i in range(1,len(m),1):
        U[0] += m[i]*(r[i][j][0] - r_mar[0]) / pow(pow(r[i][j][0] - r_mar[0], 2) + pow(r[i][j][1] - r_mar[1],2),3/2)
        U[1] += m[i]*(r[i][j][1] - r_mar[1]) / pow(pow(r[i][j][0] - r_mar[0], 2) + pow(r[i][j][1] - r_mar[1],2),3/2)
    U[0] += -G(0,j,0)
    U[1] += -G(0,j,1)
    U_mar[1] = (U[0] * np.cos(2 * np.pi * l / Z) + U[1] * np.sin(2 * np.pi * l / Z))
    U_mar[0] = -(U[0] * np.sin(2 * np.pi * l / Z) - U[1] * np.cos(2 * np.pi * l / Z))
    return [np.arctan(U_mar[1]/U_mar[0]),U[0],U[1],U_mar[0],U_mar[1]]


f = np.zeros((Z,T,D))                       #camp de vectors de força gravitatòria (cartesianes)
f_mar = np.zeros((Z,T,D))                   #camp de vectors de força gravitatòria (superfície del mar)
for i in range(Z):
    f[i%Z][0][0] = Phi((i-1)%Z,0)[1]
    f[i%Z][0][1] = Phi((i-1)%Z,0)[2]
    f_mar[i][0][0] = Phi((i-1)%Z,0)[3]
    f_mar[i][0][1] = Phi((i-1)%Z,0)[4]



for j in range(1,T-1):                          #càlcul de la posició del mar i en l'instant j+1
    for i in range(Z):
        f[i%Z][j][0] = Phi((i-1)%Z,j)[1]
        f[i%Z][j][1] = Phi((i-1)%Z,j)[2]
        f_mar[i][j][0] = Phi((i-1)%Z,j)[3]
        f_mar[i][j][1] = Phi((i-1)%Z,j)[4]

        q[i%Z][j+1] += 2*q[i%Z][j]-q[i%Z][j-1] + pow(gamma,2)*(q[(i-1)%Z][j]+q[(i+1)%Z][j]-2*q[i%Z][j])- V*dt*(q[i%Z][j]-q[i%Z][j-1])

        if (Phi(i,j)[3]>0):
            q[(i+1)%Z][j+1] += k_q*np.abs(Phi(i,j)[3])*(q[i%Z][j]/q0)*dz*pow(dt,2)
            q[(i)%Z][j+1] += -k_q*np.abs(Phi(i,j)[3])*(q[i%Z][j]/q0)*dz*pow(dt,2)
        if (Phi(i,j)[3]<0):
            q[(i-1)%Z][j+1] += k_q*np.abs(Phi(i,j)[3])*(q[i%Z][j]/q0)*dz*pow(dt,2)
            q[(i)%Z][j+1] += -k_q*np.abs(Phi(i,j)[3])*(q[i%Z][j]/q0)*dz*pow(dt,2)




#(COMPROVACIONS)
for j in range(0,T):
    S = 0                           #àrea que ocupa el mar
    CM = [0,0]                      #centre de masses del sistema
    cm = 0
    for i in range(Z):
        S += q[i][j]*dz
    for i in range(len(m)):
        cm += m[i]
        CM[0] += m[i]*r[i][j][0]
        CM[1] += m[i]*r[i][j][1]
    #print(np.asarray(CM)/cm)
    #print(S)




#ANIMACIÓ
fig = plt.figure(figsize=(15,5))
ax = fig.add_axes([1.55/30,1.25/10,8/30,8/10])          #N cossos
bx = fig.add_axes([11.25/30,1.25/10,18/30,8/10])        #mar

marc = 1.1*(np.max(r))           #marc per al plot d'òrbita

def init_func():
    ax.clear()
    bx.clear()

def update_plot(j):
    ax.clear()
    bx.clear()
    #ax.set_title("Pla d'òrbita")
    #bx.set_title("Alçada del mar")
    ax.text(0.1, 0.9, ('Temps = %.2f' % (j*dt/w)), size=10, color='black', ha='left', va='top', transform=ax.transAxes, position = (1/30,1-1/30))
    ax.text(0.1, 0.9, ('$R_M/q_0$ = %.2f' % (R[0]/q0)), size=10, color='black', ha='left', va='top', transform=ax.transAxes, position = (1/30,1/20))
    bx.text(0.1, 0.9, (r'$k_q$ = %.2f' % (k_q)), size=10, color='black', ha='left', va='top', transform=bx.transAxes, position = (1/100,1/20))
    bx.text(0.1, 0.9, (r'$w$ = %.2f' % (w)), size=10, color='black', ha='left', va='top', transform=bx.transAxes, position = (1/100,2/20))
    bx.text(0.1, 0.9, (r'$\gamma$ = %.2f' % (V)), size=10, color='black', ha='left', va='top', transform=bx.transAxes, position = (1/100,3/20))
    #ax.set(xlim=(r[0][j][0]-1.5*marc,r[0][j][0]+1.5*marc), ylim=(r[0][j][1]-1.5*marc,r[0][j][1]+1.5*marc))  #Seguiment llunyà de M
    ax.set(xlim=(r[0][j][0]-2*R[0],r[0][j][0]+2*R[0]), ylim=(r[0][j][1]-2*R[0],r[0][j][1]+2*R[0]))  #Seguiment proper de M
    #ax.set(xlim=(- marc,+ marc), ylim=(- marc,marc))                                                      #SR inercial
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    bx.set(xlim=(0,2*np.pi), ylim=(0.9*np.min(q),1.1*np.max(q)))
    bx.set_xlabel(r'$\theta$')
    bx.set_xticks(np.arange(0,2*np.pi,step=np.pi/4))
    bx.set_xticklabels(['0','π/4','π/2','3π/4','π','5π/4','3π/2', '7π/4'])
    bx.set_ylabel(r'$q$')

    print('Plotant %.2f/%.2f' % ((j*dt/w), t))

    bx.plot(np.linspace(0, 2 * np.pi, Z), q[:, j], color='b', linestyle='-', label=r'$q(\theta,t)$')
    bx.plot(np.linspace(0, 2 * np.pi, Z), q[:, 0], color='k', linestyle='--', alpha = 0.5, label=r'$q(\theta,0)$')


    Q = np.zeros((Z,2))
    Q0 = np.zeros((Z,2))
    for i in range(0,Z):
        Q[i][0] = r[0][j][0] + (R[0]+q[i][j])*np.cos(2*np.pi*(i/Z))
        Q0[i][0] = r[0][j][0] + (R[0]+q[i][0])*np.cos(2*np.pi*(i/Z))
        Q[i][1] = r[0][j][1] + (R[0]+q[i][j])*np.sin(2*np.pi*(i/Z))
        Q0[i][1] = r[0][j][1] + (R[0]+q[i][0])*np.sin(2*np.pi*(i/Z))

    ax.add_patch(Polygon(np.asarray(Q),fill=False,edgecolor='b',linewidth=2,linestyle='-'))
    ax.add_patch(Polygon(np.asarray(Q0),fill=False,edgecolor='k',alpha=0.5,linewidth=2,linestyle='--'))
    #ax.quiver(Q[:,0],Q[:,1],f[:,j,0],f[:,j,1],color='black',scale_units='xy', scale=0.0002,width=0.003)
    #bx.quiver(np.linspace(0,2*np.pi,Z),q[:,j],f_mar[:,j,0],f_mar[:,j,1],color='black',scale_units='xy',scale=0.001,label=r'$\tilde{F}(q)$',width=0.003)
    #bx.quiver(np.linspace(0,2*np.pi,Z),q[:,j],np.zeros((Z)),f_mar[:,j,1],color='red',scale_units='xy',scale=0.001,label=r'$\tilde{F}_{\perp}(q)$',width=0.001)
    #bx.quiver(np.linspace(0,2*np.pi,Z),q[:,j],f_mar[:,j,0],np.zeros((Z)),color='green',scale_units='xy',scale=0.0005,label=r'$\tilde{F}_{\parallel}(q)$',width=0.001)


    for i in range(0,len(m)):
        ax.add_patch(Circle((r[i][j][0],r[i][j][1]),R[i], fill=False, edgecolor='k', linewidth=2, linestyle='-'))

    bx.legend(loc = 'lower right')

    if j == np.max(range(0,T,step)):
        fig.savefig('marees (últim frame).png')


step = 100               #interval entre cada temps j
                                 # Atenció: Intervals petits i T grans solen donar memory leaks

anim = FuncAnimation(fig,
                     func = update_plot,
                     frames=[j for j in range(0,T,step)],
                     init_func = init_func(),
                     interval=1000
)

writergif = PillowWriter(fps=15)
anim.save('marees.gif',writer=writergif)