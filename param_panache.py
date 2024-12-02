import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

nz = 100 #nb points direction verticale
nt = 1000 #nb de pas de temps
dt = 1 #pas de temps
dz = 1 #pas vertical 
tf = nt*dt #temps final
eps = 0.005 #taux de mélange
Q = 0.05 #taux de chauffage 
alpha = 0.1 #taille du panache sur la maille


# On calcule à chaque pas de temps le b large échelle à partir des tendances 
# liées au chauffage Q et au panache paramétrisé. Chacune de ces tendances est calculée dans 
# une fonction appelée à chaque pas de temps et qui prend une entrée les grandeurs 
# grandes échelles ici b.

temps = np.linspace(1,tf+dt,nt) #temps, pas utile
b = np.zeros((nt,nz)) #b large échelle 

#Fonction de chauffage
def Qf():
    btend_Q = np.zeros(nz)
    for iz in range(nz):
        if iz>0:
            btend_Q[iz] = -Q/(nz-1)*dt
        else:
            btend_Q[iz] = Q*dt

    return btend_Q

#Fonction paramétrant le transfert lié au panache
def Panache(b):
    b_panache = np.zeros(nz) #b panache
    Gamma = np.zeros(nz) #gamma (lié à E / D)
    f = np.zeros(nz) #flux de masse entre les niveaux
    fb = np.zeros(nz) #flux de b entre les niveaux
    btend_param = np.zeros(nz) #tendance en b
    E = np.zeros(nz) #entrainement
    D = np.zeros(nz) #detrainement

    for iz in range(nz):
        if iz==0:
            #premiere couche du modèle
            f[0]=0 #flux initial nul depuis le sol

            E[iz] = alpha/2 * np.sqrt(b[0]-b[1]) * np.sqrt(dz) 
            #calcul à moitié avec le 2e niveau pour l'initialisation, on suppose f[1] = E[0] ici
            D[iz] = 0

            f[iz+1] = f[iz] + E[iz]

            b_panache[iz] = b[iz] #le b provient de l'environnement

            fb[iz+1] = f[iz+1]*(b_panache[iz]-b[iz+1])

            btend_param[iz] = -dt/dz*(fb[iz+1]-fb[iz])

           
        elif iz<nz-1:
            if f[iz] > 0:
                Gamma[iz] = alpha**2 * (b_panache[iz-1]-b[iz])/(2*f[iz])
            else:
                Gamma[iz]=0

            E[iz] = (max(Gamma[iz],0) + eps)*dz
            D[iz] = (max(-Gamma[iz],0) + eps)*dz
            
            f[iz+1] = max(f[iz] + E[iz] - D[iz],0)

            if f[iz+1]>0:
                b_panache[iz] = (f[iz] * b_panache[iz-1] + E[iz]*b[iz]) / (f[iz+1] + D[iz])
            else:
                b_panache[iz] = b_panache[iz-1]

                E[iz]=0
                D[iz] = f[iz]
            
            fb[iz+1] = f[iz+1]*(b_panache[iz]-b[iz+1])

            btend_param[iz] = -dt/dz*(fb[iz+1]-fb[iz])

        else: 
            b_panache[iz] = b_panache[iz-1]

            E[iz]=0
            D[iz] = f[iz]

            btend_param[iz] = -dt/dz*(-fb[iz])

    return btend_param
    
#profil initial de b, à importer depuis le LES de Fluids2d
b_init = np.concatenate([np.linspace(1, 0.5, nz//10, endpoint=True), np.linspace(0.5, 2, nz-nz//10, endpoint=True) ])
b[0,:] = b_init

#boucle sur les pas de temps
for it in range(nt-1):
    btend_Q = Qf() #calcul tendance en Q
    btend_param = Panache(b[it,:]) #calcul tendance panache
    b[it+1] = b[it] + btend_Q + btend_param #calcul du pas de temps suivant

# Création et écriture du fichier de sortie
filename = "output_param.nc"
ds = Dataset(filename, "w", format="NETCDF4")

ds.createDimension("t", nt)  # Dimension temporelle
ds.createDimension("z", nz)  # Dimension spatiale

b_variable = ds.createVariable("b", "f4", ("t", "z"))  # 'f4' pour float32
b_variable[:, :] = b

ds.close()


