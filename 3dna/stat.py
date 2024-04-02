import statistics as st
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import sys
from math import log2
from genetique import Genetique
from recuit import Recuit
from random import randint 
from Traj3D import Traj3D

from multiprocessing import Pool, Lock

# Read file
lineList = [line.rstrip('\n') for line in open("../data/plasmid_8k.fasta")]
# Formatting
seq = ''.join(lineList[1:])

genetique = Genetique(
    seq=seq,
    methode_utilisee="Tournoi",
    nbr_generation_max=150,
    N=50,
    probabilite_mutation_initiale=0.05,
    probabilite_mutation_finale=0.05,
    relier = 1
)

recuit = Recuit(
    seq=seq,
    k_max=1000,
    e_min=0.1,
    temp_init=300,
    refroidissement=0.99,
    dist_min=0.,
    relier=1
)

lock = Lock()

# On écrit les résultats dans un fichier
def _montecarloWriteResults(magnitude, dot, filename: str):
    with lock:
        with open(filename, 'a') as f:
            f.write(f"{magnitude},{dot}\n")

            print(f"magnitude : {magnitude} ; dot : {dot}")

# On exécute la méthode choisie
def _montecarloCompute(method, size: int):
    for _ in range(size):
        traj = Traj3D(fig = False)
        rot_table = method.executer()
        traj.compute(seq, rot_table)

        vdebut = traj.getTraj()[1] - traj.getTraj()[0]
        vmilieu = traj.getTraj()[0] - traj.getTraj()[-1]
        vfin = traj.getTraj()[-1] - traj.getTraj()[-2]

        dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
        dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

        _montecarloWriteResults(vmilieu.magnitude, min(dot1,dot2), "results.csv")

# Multiprocessing pour accélérer le calcul
def montecarloMultiProcessing(method, size: int, process: int):
    with Pool(process) as p:
        p.starmap(_montecarloCompute, [(method, size//process)]*process)
        p.close()
        p.join()

# Affiche les résultats dans un histogramme
def montecarloPlot(filename: str, nom_methode: str):
    with open(filename, 'r') as f:
        L1, L2 = [], []
        for line in f:
            magnitude, dot = line.split(',')
            L1.append(float(magnitude))
            L2.append(float(dot))
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,8))

    moy1 = st.mean(L1)
    std1 = st.stdev(L1)

    moy2 = st.mean(L2)
    std2 = st.stdev(L2)

    K = int(1 + log2(len(L1)))

    ax1.axvline(moy1, color='k', linestyle='dashed', linewidth=1, label='Moyenne')

    ax1.hist(L1, K, edgecolor='k')
    ax1.set_xlabel(f'Distance \n Moyenne = {moy1:.2f} ; Ecart type = {std1:.2f}')

    ax1.set_ylabel('Effectif')

    ax1.set_title(f"Histogramme des distances pour {len(L1)} essais, \n méthode {nom_methode}")

    ax2.axvline(moy2, color='k', linestyle='dashed', linewidth=1, label='Moyenne')

    ax2.hist(L2, K, edgecolor='k')  # Utilisation de L2 au lieu de L1
    ax2.set_xlabel(f'Minimum des similarités cos \n Moyenne = {moy2:.2f} ; Ecart type = {std2:.2f}')

    ax2.set_ylabel('Effectif')

    if nom_methode=="recuit":
        ax2.set_title(f"Histogramme des minimums des similarités cos pour {len(L2)} essais, \n méthode de recuit simulé, \n température initiale : {recuit.temp_init}, \n taux de refroidissement : {recuit.refroidissement}, kmax : {recuit.k_max}")
    else:
        ax2.set_title(f"Histogramme des minimums des similarités cos pour {len(L2)} essais, \n méthode génétique à {genetique.nbr_generation_max} gén., \n séléction par rang et p. mutation {genetique.probabilite_mutation_initiale}")

    plt.legend()
    plt.show()

def montecarlo(method, size: int) -> None:
    L1, L2 = [], [] # liste des distances, liste des PS
    
    for _ in range(size):
        traj = Traj3D(fig = False)
        rot_table = method.executer()
        traj.compute(seq, rot_table)

        vdebut = traj.getTraj()[1] - traj.getTraj()[0]
        vmilieu = traj.getTraj()[0] - traj.getTraj()[-1]
        vfin = traj.getTraj()[-1] - traj.getTraj()[-2]

        dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
        dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)
        
        L1.append(vmilieu.magnitude)
        L2.append(min(dot1,dot2))
        
    # L1 = np.random.normal(10, 3, size)
    # L2 = np.random.normal(0.1, 0.04, size)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,8))
    
    moy1 = st.mean(L1)
    var1 = st.variance(L1)
    moy2 = st.mean(L2)
    var2 = st.variance(L2)
    K = int(1 + log2(size))
    
    ax1.axvline(moy1, color='k', linestyle='dashed', linewidth=1, label='Moyenne')
    ax1.hist(L1, K, edgecolor='k')
    ax1.set_xlabel(f'Disantance \n Moyenne = {moy1:.2f} ; Variance = {var1:.2f}')
    ax1.set_ylabel('Effectif')
    if method==genetique:
        ax1.set_title(f"Histogramme des distances pour {size} essais, \n méthode géntique")
    else:
        ax1.set_title(f"Histogramme des distances pour {size} essais, \n méthode recuit simulé")

    ax2.axvline(moy2, color='k', linestyle='dashed', linewidth=1, label='Moyenne')
    ax2.hist(L2, K, edgecolor='k')  # Utilisation de L2 au lieu de L1
    ax2.set_xlabel(f'Minimum des produits scalaires \n Moyenne = {moy2:.2f} ; Variance = {var2:.2f}')
    ax2.set_ylabel('Effectif')
    if method==genetique:
        ax2.set_title(f"Histogramme des min des produits scalaires pour {size} essais, \n méthode génétique à {method.nbr_generation_max} gén., séléction par tournoi et p. mutation {method.probabilite_mutation_initiale}")
    else:
        ax2.set_title(f"Histogramme des min des produits scalaires pour {size} essais, \n méthode de recuit simulé, température initiale : {method.temp_init}, taux de refroidissement : {method.refroidissement}, kmax : {method.k_max}")
    
    plt.legend()
    plt.show()



def carte_Tk_Te(Tmin, Tmax, nT, sample_size):

    Le = []
    Lk = []
    LT = np.linspace(Tmin,Tmax,nT)
    for T in LT:
        se, sk = 0, 0
        for i in range(sample_size):
            recuit = Recuit(
            seq=seq,
            k_max=450,
            e_min=0.1,
            temp_init=T,
            refroidissement=0.99,
            dist_min=0,
            relier=1
            )

            traj = Traj3D(fig = False)
            traj.compute(seq, recuit.executer())

            vdebut = traj.getTraj()[1] - traj.getTraj()[0]
            vmilieu = traj.getTraj()[0] - traj.getTraj()[-1]
            vfin = traj.getTraj()[-1] - traj.getTraj()[-2]

            dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
            dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

            se += vmilieu.magnitude
            sk += min(dot1,dot2)
        se /= sample_size
        sk /= sample_size
        Le.append(se)
        Lk.append(sk)
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 7))
    ax1.plot(LT,Le)
    ax2.plot(LT,Lk)
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 7))
    ax1.plot(LT,Le)
    ax2.plot(LT,Lk)
    plt.tight_layout()
    plt.show()

# montecarlo(recuit, 10)

if __name__ == "__main__":
    montecarloMultiProcessing(recuit, 288, 32)
    montecarloPlot("results.csv", "recuit")