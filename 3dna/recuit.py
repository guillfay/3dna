from RotTable import RotTable
from Traj3D import Traj3D

from random import random, uniform, gauss
from math import exp

from typing import Tuple

class Recuit:
    """
    seq : La séquence de nucléotides
    k_max : Le nombre d'itération maximale
    e_min : Le seuil d'énergie
    temp_init : La température initiale
    refroidissement : 	Le coefficient de refroidissement
    dist_min : Le seuil de distance
    relier : Le nombre de nucléotides reliés entre la fin et le début de la séquence
    """
    def __init__(self,
                 seq : str,
                 k_max : int,
                 e_min : float,
                 temp_init : float,
                 refroidissement : float,
                 dist_min : float,
                 relier : int = 0):
        self.seq = seq
        self.k_max = k_max
        self.e_min = e_min
        self.temp_init = temp_init
        self.refroidissement = refroidissement
        self.dist_min = dist_min
        self.relier = relier

        # On garde une copie de la version originale de RotTable
        self.rot_table = RotTable()

    def energie(self, s : RotTable) -> Tuple[float, float]:
        traj = Traj3D(fig = False)
        traj.compute(self.seq, s)

        if self.relier > 0:
            sommeDist = 0

            # On somme la distance entre les nucléotides à relier
            for i in range(self.relier):
                debut = traj.getTraj()[i]
                fin = traj.getTraj()[-self.relier + i]

                sommeDist += (debut - fin).magnitude
            
            # vdebut est le vecteur entre les deux premier nucléotides
            # vmilieu est le vecteur entre le dernier et premier nucléotide
            # vfin est le vecteur entre l'avant dernier et le dernier nucléotide
            vdebut = traj.getTraj()[1] - traj.getTraj()[0]
            vmilieu = traj.getTraj()[0] - traj.getTraj()[-self.relier-1]
            vfin = traj.getTraj()[-self.relier-1] - traj.getTraj()[-self.relier-2]
                
            dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
            dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

            dist = (traj.getTraj()[0] - traj.getTraj()[-self.relier - 1]).magnitude
            
            # L'énergie est la somme entre la moyenne des distances entre les nucléotides à relier, et un facteur possédant une valeur élevée lorsque les produits scalaires entre les trois vecteurs s'écarte de 1.
            return (sommeDist / self.relier) + 6. * exp(-dist / 50) * (1. - exp(-0.8*((dot1 - 1) ** 2 + (dot2 - 1) ** 2))), dist

        debut = traj.getTraj()[0]
        fin = traj.getTraj()[-1]

        dist = (debut - fin).magnitude

        # vdebut est le vecteur entre les deux premier nucléotides
        # vmilieu est le vecteur entre le dernier et premier nucléotide
        # vfin est le vecteur entre l'avant dernier et le dernier nucléotide
        vdebut = traj.getTraj()[1] - traj.getTraj()[0]
        vmilieu = traj.getTraj()[0] - traj.getTraj()[-1]
        vfin = traj.getTraj()[-1] - traj.getTraj()[-2]

        dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
        dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

        # L'énergie est la somme entre différence entre la distance et 3, et un facteur possédant une valeur élevée lorsque les produits scalaires entre les trois vecteurs s'écarte de 1.
        return abs(dist - 3) + 5. * exp(-dist / 50) * (1. - exp(-0.8*((dot1 - 1) ** 2 + (dot2 - 1) ** 2))), dist

    def voisin(self, s : RotTable, t : float) -> RotTable:
        # On copie le RotTable actuel
        sn = RotTable(rot_table = s.rot_table)

        # On divise la variance de la gaussienne par 3, et elle décroît avec la température qui descend.
        facteur_std = (t / self.temp_init) / 3.

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
        # Pour chaque nucléotide, on modifie la valeur en utilisant une distribution gaussienne qui a pour espérance la valeur antérieure, et pour écart type le produit entre le facteur précédent et la limite.
        # On copie ces valeurs à leur complémentaire s'il existe
        # for key in sn.rot_table:
        for (key, key2) in PAIRES:
            sn.rot_table[key][0] = gauss(sn.rot_table[key][0], facteur_std * sn.rot_table[key][3])
            sn.rot_table[key][1] = gauss(sn.rot_table[key][1], facteur_std * sn.rot_table[key][4])

            # On borne les différentes valeurs
            sn.rot_table[key][0] = max(
                self.rot_table.rot_table[key][0] - self.rot_table.rot_table[key][3],
                min(
                    self.rot_table.rot_table[key][0] + self.rot_table.rot_table[key][3],
                    sn.rot_table[key][0]
                )
            )

            sn.rot_table[key][1] = max(
                self.rot_table.rot_table[key][1] - self.rot_table.rot_table[key][4],
                min(
                    self.rot_table.rot_table[key][1] + self.rot_table.rot_table[key][4],
                    sn.rot_table[key][1]
                )
            )

            if key2 is not None:
                sn.rot_table[key2][0] = sn.rot_table[key][0]
                sn.rot_table[key2][1] = sn.rot_table[key][1]
        
        return sn

    def P(self, e_delta: float, t: float) -> float:
        if e_delta < 0:
            return 1.

        return exp(-e_delta / t)
    
    # Algorithme de recuite simulée
    def executer(self) -> RotTable:
        if self.relier > 0:
            print("Energie : relier")

            self.seq += self.seq[:self.relier]

        s = RotTable()
        e, distance = self.energie(s)
        k = 0
        temp = self.temp_init

        while k < self.k_max and e > self.e_min and distance > self.dist_min:
            sn = self.voisin(s, temp)
            en, distancen = self.energie(sn)

            if en < e or random() < self.P(en - e, temp):
                s = RotTable(rot_table=sn.rot_table)
                e = en
                distance = distancen

            print(f"{k}\tdist : {distance}\te : {e}\ten : {en}\ttemp : {temp}")

            k += 1
            temp *= self.refroidissement

        if self.relier > 0:
            self.seq = self.seq[:-self.relier]

        return s