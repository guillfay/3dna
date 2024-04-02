from random import random
import random
import sys
import matplotlib.pyplot as plt
import numpy as np
from recuit import  Recuit
from genetique import Genetique
from RotTable import RotTable
from Traj3D import Traj3D

from multiprocessing import Pool, Lock

# Réécriture de la fonction exécuter de Recuit()
# Renvoie en plus la liste des énergies pour faire un plot
def executer_plot(recuit):
    if recuit.relier > 0:
        print("Energie : relier")

        recuit.seq += recuit.seq[:recuit.relier]

        s = RotTable()
        e, distance = recuit.energie(s)
        k = 0
        temp = recuit.temp_init
        e_liste=[]
        while k < recuit.k_max and e > recuit.e_min and distance > recuit.dist_min:
            sn = recuit.voisin(s, temp)
            en, distancen = recuit.energie(sn)

            if en < e or random() < recuit.P(en - e, temp):
                s = RotTable(rot_table=sn.rot_table)
                e = en
                distance = distancen

            print(f"{k}\tdist : {distance}\te : {e}\ten : {en}\ttemp : {temp}")
            e_liste.append(e)
            k += 1
            temp *= recuit.refroidissement

        if recuit.relier > 0:
            recuit.seq = recuit.seq[:-recuit.relier]

        return s,e_liste

# Réécriture de la fonction exécuter de Genetique()
# Renvoie en plus la liste des énergies pour faire un plot    
def executer_plot_genetique(genetique):
    if genetique.relier > 0:
        print("Energie : relier")

        genetique.seq += genetique.seq[:genetique.relier]

    s = RotTable()
    e_liste=[]
    generation = 0
    genetique.probabilite_mutation = genetique.probabilite_mutation_initiale
    population = genetique.generation_population()

    # Selection de la méthode utilisée pour la selection
    if genetique.methode_utilisee == "elitisme":
        methode_selection = genetique.selection_elitisme
        print("\nMethode utilisée : Elitisme\n")
    elif genetique.methode_utilisee == "roulette":
        methode_selection = genetique.selection_roulette
        print("\nMethode utilisée : Roulette\n")
    elif genetique.methode_utilisee == "rang":
        methode_selection = genetique.selection_rang
        print("\nMethode utilisée : Rang\n")
    elif genetique.methode_utilisee == "tournoi":
        methode_selection = genetique.selection_tournoi
        print("\nMethode utilisée : Tournoi\n")
    else:  # Par défaut
        methode_selection = genetique.selection_elitisme
        print("\nMethode utilisée par défaut : Elitisme\n")

    # Lancement de la boucle génétique
    while generation < genetique.nbr_generation_max:
        population_selectionnes = methode_selection(population)

        # Mutation du/des max, c'est un agent "crash-test" explorateur
        genetique.mutation(population_selectionnes[-1], 0.9, 0.04)
        genetique.mutation(population_selectionnes[-2], 0.06, 0.5)

        if generation % 50 == 0:
            print(
                f"{generation}\te_min : {round(genetique.energie(population_selectionnes[0]),5)}\te_moy : {round(sum([genetique.energie(population_selectionnes[i]) for i in range(len(population_selectionnes))])/len(population_selectionnes), 5)}\te_max : {round(genetique.energie(population_selectionnes[-1]), 5)}\tproba_mutation : {round(genetique.probabilite_mutation, 5)}")
        e_liste.append(round(genetique.energie(population_selectionnes[0]),5))
        # Croisement des individus
        random.shuffle(population_selectionnes)
        population_enfant = []
        for i in range(0, len(population_selectionnes)-1, 2):
            population_enfant += genetique.croisement_N_points(
                population_selectionnes[i], population_selectionnes[i+1], 13)

        # Reconstitution population
        population = population_selectionnes + population_enfant

        population = sorted(population, key=lambda x: genetique.energie(x))

        # Mutations
        for i in range(1, len(population)):
            genetique.mutation(
                population[i], (1 - ((generation + 1) / genetique.nbr_generation_max)) * 0.33, genetique.probabilite_mutation)

        # Modification de la probabilité
        genetique.probabilite_mutation = genetique.probabilite_mutation_finale*generation / \
            genetique.nbr_generation_max + \
            genetique.probabilite_mutation_initiale*(genetique.nbr_generation_max-generation) / \
            genetique.nbr_generation_max

        # Génération suivante
        generation += 1

    if genetique.relier > 0:
        genetique.seq = genetique.seq[:-genetique.relier]

    return min(population, key=lambda x: genetique.energie(x)), e_liste

# Recuit : Energie en fonction du nombre d'itérations pour plasmid_8k
def recuit_e_f_k_8k(T,nb_sim):
    e_global=[]
    for k in range(nb_sim):
        traj = Traj3D(fig=False)

        seq=''.join([line.rstrip('\n') for line in open('../data/plasmid_8k.fasta')][1:])

        recuit = Recuit(
            seq = seq,
            k_max = 700,
            e_min = 0.1,
            temp_init = T,
            refroidissement = 0.99,
            dist_min=0.,
            relier=1
            )
     
        rot_table,e = executer_plot(recuit)
        e_global.append(e)
    fig=plt.figure()
    ax=plt.axes()
    print(e_global)
    for e in e_global:
        plt.plot([k for k in range(len(e))],e, label=str(round(e[-1],2)))
    plt.title(f"{str(nb_sim)} simulations à T_init={str(T)}")
    plt.xlabel("Nombre d'itérations")
    plt.ylabel("Energie")
    plt.legend()
    plt.show()

# Recuit : Energie en fonction du nombre d'itérations pour plasmid_180k
def recuit_e_f_k_180k(T,nb_sim):
    e_global=[]
    for k in range(nb_sim):
        traj = Traj3D(fig=False)

        seq=''.join([line.rstrip('\n') for line in open(r'..\data\plasmid_180k.fasta')][1:])

        recuit = Recuit(
            seq = seq,
            k_max = 700,
            e_min = 0.1,
            temp_init = T,
            refroidissement = 0.99,
            dist_min=0.,
            relier=1)
            # seq = seq,
            # k_max = 500,
            # e_max = 1,
            # temp_init = 300,
            # refroidissement = 0.98 ==> e=0.28
        

        rot_table,e = executer_plot(recuit)
        e_global.append(e)
    fig=plt.figure()
    ax=plt.axes()
    print(e_global)
    for e in e_global:
        plt.plot([k for k in range(len(e))],e, label=str(round(e[-1],2)))
    plt.title(f"{str(nb_sim)} simulations à T_init={str(T)}")
    plt.xlabel("Nombre d'itérations")
    plt.ylabel("Energie")
    plt.legend()
    plt.show()

# Recuit : Energie en fonction du nom pour plasmid_8k
def recuit_e_f_T(nb_sim):
    e_global=[]
    fig=plt.figure()
    ax=plt.axes()
    for t in range(100,1100,100):
        e_local=[]
        for k in range(nb_sim):
            traj = Traj3D(fig=False)

            seq=''.join([line.rstrip('\n') for line in open('../data/plasmid_8k.fasta')][1:])

            recuit = Recuit(
                seq = seq,
                k_max = 500,
                e_min = 0.1,
                temp_init = t,
                refroidissement = 0.99,
                dist_min=0.,
                relier=1)
                # seq = seq,
                # k_max = 500,
                # e_max = 1,
                # temp_init = 300,
                # refroidissement = 0.98 ==> e=0.28
            

            rot_table,e = executer_plot(recuit)
            e_local.append(e[-1])
        plt.scatter([t]*len(e_local),e_local)
        e_global.append(e_local)
    
    print(e_global)
    plt.title("Répartition Energie VS Temp_init")
    plt.xlabel('Température initiale')
    plt.ylabel('Energie')
    plt.legend()
    plt.show()
    
lock = Lock()

def compute_write_results(n, e):
    with lock:
        with open("results.csv", 'a') as f:
            f.write(f"{n},{e}\n")
            print(f"n : {n} ; e : {e}")

def compute_genetique_e_f_N(i):
    traj = Traj3D(fig=False)

    seq=''.join([line.rstrip('\n') for line in open('../data/plasmid_8k.fasta')][1:])

    for n in range(20, 121, 20):
        genetique = Genetique(
            seq=seq,
            methode_utilisee='elitisme',
            nbr_generation_max=150,
            N=n,
            probabilite_mutation_initiale=0.05,
            probabilite_mutation_finale=0.05,
            relier = 1
        )

        rot_table,e = executer_plot_genetique(genetique)
        
        compute_write_results(n, e[-1])

def multiprocessing_genetique_e_f_N(thread):
    with Pool(thread) as p:
        p.map(compute_genetique_e_f_N, [0 for _ in range(thread)])

def genetique_e_f_N_plot(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        n = [int(line.split(',')[0]) for line in lines]
        e = [float(line.split(',')[1]) for line in lines]
    
    fig=plt.figure()
    ax=plt.axes()

    last_n = -1
    e_local = []
    for i in range(len(n)):
        if n[i] != last_n:
            if len(e_local) > 0:
                plt.scatter([last_n]*len(e_local), e_local)

                e_np = np.array(e_local)

                print(f"n = {last_n} Moyenne : {e_np.mean()} ; Ecart-type : {e_np.std()}")

            e_local = []
            last_n = n[i]
        
        e_local.append(e[i])
    
    plt.scatter([last_n]*len(e_local), e_local)

    e_np = np.array(e_local)

    print(f"n = {last_n} Moyenne : {e_np.mean()} ; Ecart-type : {e_np.std()}")

    plt.title("Répartition Energie VS Nombre d'individus")
    plt.xlabel('Nombre individus')
    plt.ylabel('Energie')
    plt.ylim(0,20)
    # plt.legend()
    plt.show()

def genetique_e_f_N(nb_sim):
    e_global=[]
    fig=plt.figure()
    ax=plt.axes()
    for n in range(20,100,20):
        e_local=[]
        for k in range(nb_sim):
            traj = Traj3D(fig=False)

            seq=''.join([line.rstrip('\n') for line in open('../data/plasmid_8k.fasta')][1:])

            genetique = Genetique(
                seq=seq,
                methode_utilisee='elitisme',
                nbr_generation_max=150,
                N=n,
                probabilite_mutation_initiale=0.05,
                probabilite_mutation_finale=0.05,
                relier = 1
                # seq = seq,
                # k_max = 500,
                # e_min = 0.1,
                # temp_init = t,
                # refroidissement = 0.99,
                # dist_min=0.,
                # relier=1
                )
                # seq = seq,
                # k_max = 500,
                # e_max = 1,
                # temp_init = 300,
                # refroidissement = 0.98 ==> e=0.28
            

            rot_table,e = executer_plot_genetique(genetique)
            e_local.append(e[-1])
        plt.scatter([n]*len(e_local),e_local)
        e_global.append(e_local)
    
    print(e_global)
    plt.title("Répartition Energie VS Nombre d'individus")
    plt.xlabel('Nombre individus')
    plt.ylabel('Energie')
    plt.ylim(0,50)
    plt.legend()
    plt.show()

# genetique_e_f_N(2)
multiprocessing_genetique_e_f_N(32)
genetique_e_f_N_plot("results.csv")