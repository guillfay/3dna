from RotTable import RotTable
from Traj3D import Traj3D

import random
from math import exp
import copy

from typing import List, Callable, Tuple

rot_table = RotTable()

PAIRES = [
    ("AA", "TT"),
    ("AC", "GT"),
    ("AG", "CT"),
    ("AT", None),
    ("CA", "TG"),
    ("CC", "GG"),
    ("CG", None),
    ("GA", "TC"),
    ("GC", None),
    ("TA", None)
]

"""
seq : La séquence de nucléotides
methode_utilisee : Le nom de la stratégie de sélection
nbr_generation_max : Le nombre de génération maximal
N : La taille de la population
probabilite_mutation_initiale : La probabilité de mutation initiale
probabilite_mutation_finale : La probabilité de mutation finale
relier : Le nombre de nucléotides reliés entre la fin et le début de la séquence
"""
class Genetique:
    def __init__(self,
                 seq: str,
                 methode_utilisee: str,
                 nbr_generation_max: int,
                 N: int,
                 probabilite_mutation_initiale: float,
                 probabilite_mutation_finale : float,
                 relier : int):
        self.seq = seq
        self.methode_utilisee = methode_utilisee
        self.nbr_generation_max = nbr_generation_max
        # Nombre individus
        self.N = N
        # Probaiblité de mutation génération 1
        self.probabilite_mutation_initiale = probabilite_mutation_initiale
        # Probaiblité de mutation génération finale
        self.probabilite_mutation_finale = probabilite_mutation_finale

        self.relier = relier

        self.energie_memo = {}

    # On utilise la même énergie que dans l'algorithme du recuit simulé, mais on utilise la mémoisatiion afin de mettre les valeurs en cache.
    def energie(self, s: RotTable) -> float:
        if s in self.energie_memo:
            return self.energie_memo[s]

        traj = Traj3D(fig = False)
        traj.compute(self.seq, s)

        if self.relier > 0:
            sommeDist = 0

            for i in range(self.relier):
                debut = traj.getTraj()[i]
                fin = traj.getTraj()[-self.relier + i]

                sommeDist += (debut - fin).magnitude
            
            vdebut = traj.getTraj()[1] - traj.getTraj()[0]
            vmilieu = traj.getTraj()[0] - traj.getTraj()[-self.relier-1]
            vfin = traj.getTraj()[-self.relier-1] - traj.getTraj()[-self.relier-2]
                
            dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
            dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

            dist = (traj.getTraj()[0] - traj.getTraj()[-self.relier - 1]).magnitude
            
            out = (sommeDist / self.relier) + 6. * exp(-dist / 50) * (1. - exp(-0.8*((dot1 - 1) ** 2 + (dot2 - 1) ** 2)))

            self.energie_memo[s] = out

            return out

        debut = traj.getTraj()[0]
        fin = traj.getTraj()[-1]

        dist = (debut - fin).magnitude

        vdebut = traj.getTraj()[1] - traj.getTraj()[0]
        vmilieu = traj.getTraj()[0] - traj.getTraj()[-1]
        vfin = traj.getTraj()[-1] - traj.getTraj()[-2]

        dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
        dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

        out = abs(dist - 3) + 5. * exp(-dist / 50) * (1. - exp(-0.8*((dot1 - 1) ** 2 + (dot2 - 1) ** 2)))

        self.energie_memo[s] = out

        return out

    # On génére les individus selon une loi gaussienne
    def generation_individu(self):
        table_individu = RotTable()

        for (cle, cle2) in PAIRES:
            for i in range(2):
                valeur_min = table_individu.rot_table[cle][i] - \
                    table_individu.rot_table[cle][i+3]
                valeur_max = table_individu.rot_table[cle][i] + \
                    table_individu.rot_table[cle][i+3]
                # Mise à jour de la valeur dans la table
                # table_individu.rot_table[cle][i] = random.uniform(
                #     valeur_min, valeur_max)

                table_individu.rot_table[cle][i] = min(max(random.gauss(
                    table_individu.rot_table[cle][i], table_individu.rot_table[cle][i+3] / 3
                ), valeur_min), valeur_max)

                if cle2 is not None:
                    table_individu.rot_table[cle2][i] = table_individu.rot_table[cle][i]

        return table_individu

    def generation_population(self):
        liste_individu = []

        # On génére N individus
        for i in range(self.N):
            liste_individu.append(self.generation_individu())
        return liste_individu

    def selection_elitisme(self, population: List[RotTable]) -> List[RotTable]:
        population.sort(key=lambda x: self.energie(x))

        # On ne garde que la première moitié
        return population[: len(population) // 2]

    def selection_roulette(self, population: List[RotTable]) -> List[RotTable]:
        # Inverser les scores
        # Inverse car notre score attribue de hautes valeurs aux mauvais individus
        scores_inverses = [1/self.energie(individu)for individu in population]
        somme_scores = sum(scores_inverses)
        probabilites = [score_individu /
                        somme_scores for score_individu in scores_inverses]

        selections = []
        for _ in range(len(population) // 2):
            choix = random.choices(population, weights=probabilites, k=1)[0]
            selections.append(choix)

        return sorted(selections, key=lambda x: self.energie(x))

    def selection_rang(self, population: List[RotTable]) -> List[RotTable]:
        # reversed car notre score attribue de hautes valeurs aux mauvais individus
        population_triee = sorted(
            population, key=lambda x: self.energie(x), reverse=True)
        probabilites = [(i+1)/len(population_triee)
                        for i in range(len(population_triee))]
        selections = []
        for _ in range(len(population) // 2):
            choix = random.choices(
                population_triee, weights=probabilites, k=1)[0]
            selections.append(choix)

        return sorted(selections, key=lambda x: self.energie(x))

    def selection_tournoi(self, population: List[RotTable]) -> List[RotTable]:
        population_gagnante = []

        p = 0.01

        while len(population_gagnante) < len(population) // 2 - 1 :
            individus = random.sample(population, k=2)
            scores = [self.energie(individu) for individu in individus]

            if scores[0] < scores[1]:
                score_gagnant, score_perdant = scores[0], scores[1]
                gagnant, perdant = individus[0], individus[1]
            else:
                score_gagnant, score_perdant = scores[1], scores[0]
                gagnant, perdant = individus[1], individus[0]

            if random.random() < p or self.methode_utilisee=='test_tournoi_p':
                population_gagnante.append((score_perdant, perdant))
            else:
                population_gagnante.append((score_gagnant, gagnant))

        # Selection du dernier individu

        return [min(population, key=lambda x: self.energie(x))] + [x[1] for x in sorted(population_gagnante, key=lambda x: x[0])]

    def croisement(self, individu1: RotTable, individu2: RotTable) -> [RotTable, RotTable]:
        enfant_1 = RotTable(rot_table=individu1.rot_table)
        enfant_2 = RotTable(rot_table=individu2.rot_table)

        for i in range(random.randint(0, 16)):
            k = list(individu2.rot_table.keys())[i]

            enfant_1.rot_table[k] = individu2.rot_table[k].copy()
            enfant_2.rot_table[k] = individu1.rot_table[k].copy()

        return [enfant_1, enfant_2]

    # On prend N points de coupes aléatoires. A chaque fois qu'on passe un point de coupe, on change de parent.
    def croisement_N_points(self, individu1: RotTable, individu2: RotTable, N: int) -> [RotTable, RotTable]:
        if N > len(individu1.rot_table.keys()) - 1:
            raise Exception("Nombre de coupe trop élevée")

        enfant_1 = RotTable(rot_table=individu1.rot_table)
        enfant_2 = RotTable(rot_table=individu2.rot_table)

        cles = list(individu2.rot_table.keys())

        idx = list(range(len(cles)))
        random.shuffle(idx)

        points_de_coupe = list(sorted(idx[:N]))

        for j in range(len(points_de_coupe)):
            inverser = True

            dernier_indice = points_de_coupe[j + 1] \
                if j < len(points_de_coupe) - 1 else len(cles) - 1

            if inverser:
                for i in range(points_de_coupe[j], dernier_indice + 1):
                    enfant_1.rot_table[cles[i]
                                       ] = individu2.rot_table[cles[i]].copy()
                    enfant_2.rot_table[cles[i]
                                       ] = individu1.rot_table[cles[i]].copy()

            inverser = not inverser

        return [enfant_1, enfant_2]

    def mutation(self, individu: RotTable, stdFact: float = 1., probabilite_mutation: float = 1.):
        for (key, key2) in PAIRES:
            if random.random() < probabilite_mutation:
                for i in range(2):
                    # individu.rot_table[key][1] = random.uniform(
                    #     rot_table.rot_table[key][1] - rot_table.rot_table[key][4],
                    #     rot_table.rot_table[key][1] + rot_table.rot_table[key][4]
                    # )
                    valeur_min = rot_table.rot_table[key][i] - rot_table.rot_table[key][i + 3]
                    valeur_max = rot_table.rot_table[key][i] + rot_table.rot_table[key][i + 3]
                
                    individu.rot_table[key][i] = min(max(random.gauss(
                        individu.rot_table[key][i], rot_table.rot_table[key][i + 3] * stdFact
                    ), valeur_min), valeur_max)

                if key2 is not None:
                    individu.rot_table[key2][0] = rot_table.rot_table[key][0]
                    individu.rot_table[key2][1] = rot_table.rot_table[key][1]

    def executer(self) -> RotTable:
        if self.relier > 0:
            print("Energie : relier")

            self.seq += self.seq[:self.relier]

        s = RotTable()
        generation = 0
        self.probabilite_mutation = self.probabilite_mutation_initiale
        population = self.generation_population()

        # Selection de la méthode utilisée pour la selection
        if self.methode_utilisee == "elitisme":
            methode_selection = self.selection_elitisme
            print("\nMethode utilisée : Elitisme\n")
        elif self.methode_utilisee == "roulette":
            methode_selection = self.selection_roulette
            print("\nMethode utilisée : Roulette\n")
        elif self.methode_utilisee == "rang":
            methode_selection = self.selection_rang
            print("\nMethode utilisée : Rang\n")
        elif self.methode_utilisee == "tournoi":
            methode_selection = self.selection_tournoi
            print("\nMethode utilisée : Tournoi\n")
        else:  # Par défaut
            methode_selection = self.selection_elitisme
            print("\nMethode utilisée par défaut : Elitisme\n")

        # Lancement de la boucle génétique
        while generation < self.nbr_generation_max:
            population_selectionnes = methode_selection(population)

            # Mutation du/des max, c'est un agent "crash-test" explorateur
            self.mutation(population_selectionnes[-1], 0.9, 0.04)
            self.mutation(population_selectionnes[-2], 0.06, 0.5)

            print(
                f"{generation}\te_min : {round(self.energie(population_selectionnes[0]),5)}\te_moy : {round(sum([self.energie(population_selectionnes[i]) for i in range(len(population_selectionnes))])/len(population_selectionnes), 5)}\te_max : {round(self.energie(population_selectionnes[-1]), 5)}\tproba_mutation : {round(self.probabilite_mutation, 5)}")

            # Croisement des individus
            random.shuffle(population_selectionnes)
            population_enfant = []
            for i in range(0, len(population_selectionnes)-1, 2):
                population_enfant += self.croisement_N_points(
                    population_selectionnes[i], population_selectionnes[i+1], 13)

            # Reconstitution population
            population = population_selectionnes + population_enfant

            population = sorted(population, key=lambda x: self.energie(x))

            # Mutations
            for i in range(1, len(population)):
                self.mutation(
                    population[i], (1 - ((generation + 1) / self.nbr_generation_max)) * 0.33, self.probabilite_mutation)

            # Modification de la probabilité
            self.probabilite_mutation = self.probabilite_mutation_finale*generation / \
                self.nbr_generation_max + \
                self.probabilite_mutation_initiale*(self.nbr_generation_max-generation) / \
                self.nbr_generation_max

            # Génération suivante
            generation += 1
        
        if self.relier > 0:
            self.seq = self.seq[:-self.relier]

        return min(population, key=lambda x: self.energie(x))
