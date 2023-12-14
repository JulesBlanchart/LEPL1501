
import math
import matplotlib.pyplot as plt
import numpy as np

### Constantes
g = 9.81 	# gravitation [m/s**2]
rho= 1000 	#masse volumique de l'eau [kg/m**3]
D = 0.70 	#coefficient déstabilisateur/ perte d'énergie
rad= 180/math.pi #convertir les radians en degrés

### Variables
#barge:
h1= 0.06 #hauteur [m]
m1= 1.18 #masse [kg]
L1= 0.6 #longueur [m]
Cg1 = (0,0) #centre gravité, hypothèse: la distance (d**2) entre le centre de rotation et le centre de masse 1 est négligeable

#mat:
h2= 0.43 #hauteur [m]
m2= 1.12 #masse [kg]
L2 = 0.05 #longueur [m]
Cg2=(0, h1 + h2/2) 	#centre de gravité (x,y)


#bras:
h3= 0.050 #hauteur [m]
m3= 0.95 #masse [kg]
L3= 0.95 #longueur [m]
x3=  0.21 #longueur du bras à l'arriere de la grue
Cg3=(L3/2 - x3, h1+ h2 +h3/2) #centre de gravité (x,y)

#contre-poids:
m4= 0.75 #masse [m]
Cg4= (-x3, h1+ h2 +h3/2)

#déplacement:
dx= 0.71 #déplacement en x [m]
dy= 0.54 #hauteur de la masse déplacée [m]
d= math.sqrt(dx**2 + dy**2)#calcul de la distance entre centre de rotation et la position du déplacement
md= 0.55 # masse grappin [kg]

###calculs
m_tot = m1 + m2 + m3 +m4 +md #masse totale
h_im = m_tot/(rho*(L1**2)) #hauteur immergée
angle_sub = math.atan(2*(h1-h_im)/L1) #angle d'inclinaison critique, calcul de submersion
angle_soul = math.atan(2*h_im/L1) #angle de soulèvement
I= 1.25165 #inertie totale


### Paramètres de la simulation

step = 0.001     # pas (dt) [s]
end = 15.0       # durée [s]
alpha_0 = 0.037  # angle initiale [rad] = 2,12 degrés
#0.017
w_0 = 0.0        # vitesse angulaire [m/s]


t = np.arange(0, end, step)
alpha = np.empty_like(t)
w = np.empty_like(t)
a = np.empty_like(t)
E_k = np.empty_like(t) #energie cinétique
E_p = np.empty_like(t) #energie poussée
E_t = np.empty_like(t) #energie totale

###fonctions

def simulation(d,m):
    """
    pre: -
    post: exécute une simulation jusqu'à t=end par pas de dt=step.
          Remplit les listes x, v, a des positions, vitesses et accélérations.
    """
    #conditions initiales
    alpha[0] = alpha_0
    w[0]= w_0
    E_k[0] = 0
    E_p[0] = 0
    E_t[0] = 0
    for i in range(len(t)-1):
        
        dt = step
        
        #calculs des centres de gravité et poussée
            #centre de gravité
        Cg_x = ((m1 * Cg1[0] + m2 * Cg2[0] + m3 * Cg3[0] + Cg4[0]*m4 ) / (m_tot))
       
            #centre de poussée
        Cp_x = (L1**2 *math.tan(alpha[i])/(12*h_im))
        
        #calcul des couples
        
        Cr= (m_tot) * g * abs(Cp_x - Cg_x) #couple de redressement 
        Ca= -m*g*d #couple déstabilisateur
        Cd= -D * w[i] #couple d'amortissement(a<)
        
        somme_couple =  Ca + Cr + Cd
        
        #calculs accélération, vitesse angulaire, angle
        a[i]=somme_couple / I #accélération
        w[i+1]=w[i]+a[i]*dt #vitesse angulaire
        alpha[i+1]= alpha[i] + w[i]*dt #angle
        a[i+1]=a[i]
        
        # Calculs des énergies
        
        E_k[i+1]= (I*m_tot*w[i]**2)/2
        E_p[i+1] = -Ca * alpha[i]
        E_t[i+1] = E_k[i+1] +E_p[i+1]
        

###graphiques:
        
###graphique représentant l'angle, la vitesse angulaire et l'accélération en fonction du temps [s]
def graphiques():
   
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
 
    #(mesure_t, mesure_x, mesure_v) = np.loadtxt("mesure(6).txt").T
 
    # 1er subplot
    ax1.plot(t,alpha, label="alpha")
    #ax1.plot(mesure_t, mesure_x, label = "angle mesure")
    ax1.plot([0,end], [angle_sub,angle_sub], '--r', label='submersion')    # submersion (positif)
    ax1.plot([0,end], [-angle_sub,-angle_sub], '--r')                      # submersion (négatif)
    ax1.plot([0,end], [angle_soul,angle_soul], '--b', label='soulèvement') # soulèvement (positif)
    ax1.plot([0,end], [-angle_soul,-angle_soul], '--b')                    # soulèvement (négatif)
    ax1.set_ylabel('angle [rad]')
    ax1.legend()
 
    # 2e subplot
    ax2.plot(t,w, label="vitesse angulaire")
    #ax2.plot(mesure_t, mesure_v,label = "vitesse mesure")
    ax2.set_ylabel('vitesse [rad/s]')
    ax2.legend()
 
    # 3e subplot
    ax3.plot(t, a, label='acceleration')
    ax3.set_xlabel('temps [s]')
    ax3.set_ylabel('acc [rad/s²]')
    ax3.legend()
 
    fig.suptitle("Oscillations, vitesse, accélération") # pour mettre un grand titre tout au dessus
    plt.show()
    

        
### Diagramme de phase : graphique représentant la vitesse angulaire [rad/s] en fonction de l'angle [rad]
 
def diagramme_de_phase():
    # Paramètres généraux du graphique
    
    simulation(dx,md)
    plt.title("Diagramme de phase")
    plt.plot(alpha*rad, w*rad,linewidth=0.75)
    plt.xlabel("Angle (°)")
    plt.ylabel("Vitesse angulaire (°/s)")
    plt.show()

### Graphique des Energies: graphique représentant l'énergie cinétique, de poussée et totale
def graphiques_energie():
    
    simulation(dx,md)
    plt.figure("Energies")
    plt.plot(t,E_k, "g", label="Energie cinétique")
    plt.plot(t, E_p, "r", label="Energie de poussée")
    plt.plot(t, E_t, "m", label = "Energie totale")
    plt.xlabel("temps [s]")
    plt.ylabel("Energie [J]")
    plt.suptitle("Energies")
    plt.legend()
    plt.show()
    

simulation(dx,md)
graphiques()
diagramme_de_phase()
graphiques_energie()
