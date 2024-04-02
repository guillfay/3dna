from RotTable import RotTable
from Traj3D import Traj3D
from recuit import Recuit
from genetique import Genetique

import json
import argparse

# On récupère les arguments passés en ligne de commande
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="input filename of DNA sequence")

# On ajoute un argument pour choisir l'algorithme à utiliser
parser.add_argument('mode')

# On ajoute un argument pour choisir la méthode utilisée pour la sélection de l'algorithme génétique
parser.add_argument('methode_utilisee', nargs='?', default=None,
                    help='Description de la méthode utilisée pour la sélection de l"algorithme génétique')
parser.parse_args()
args = parser.parse_args()


def main():

    traj = Traj3D()

    # Read file
    lineList = [line.rstrip('\n') for line in open(args.filename)]
    # Formatting
    seq = ''.join(lineList[1:])

    # On exécute l'algorithme choisi
    if args.mode == "recuit":
        recuit = Recuit(
            seq=seq,
            k_max=1000,
            e_min=0.1,
            temp_init=300,
            refroidissement=0.99,
            dist_min=0.,
            relier=1
        )

        rot_table = recuit.executer()

    elif args.mode == "genetique":
        genetique = Genetique(
            seq=seq,
            methode_utilisee=args.methode_utilisee,
            nbr_generation_max=150,
            N=50,
            probabilite_mutation_initiale=0.05,
            probabilite_mutation_finale=0.05,
            relier = 1
        )

        rot_table = genetique.executer()

    # On affiche le nouveau json de RotTable
    print("\nTable des rotations 3D :\n", json.dumps(rot_table.rot_table))

    # On reprend le json original pour comparer les valeurs
    rot_table_original = RotTable()

    print("\nLa différence entre les valeurs de Twist et Wedge et leur seuil maximal pour chaque paire de nucléotides (entre parenthèses) :")
    # On affiche les différences des valeurs
    for key in rot_table_original.rot_table:
        delta1 = abs(
            rot_table_original.rot_table[key][0] - rot_table.rot_table[key][0])
        delta2 = abs(
            rot_table_original.rot_table[key][1] - rot_table.rot_table[key][1])

        print(
            f"{key}\t{delta1:.4f} ({rot_table_original.rot_table[key][3]})\t{delta2:.4f} ({rot_table_original.rot_table[key][4]})")

    # On calcule la trajectoire
    traj.compute(seq, rot_table)

    # vdebut est le vecteur entre les deux premier nucléotides
    # vmilieu est le vecteur entre le dernier et premier nucléotide
    # vfin est le vecteur entre l'avant dernier et le dernier nucléotide
    vdebut = traj.getTraj()[1] - traj.getTraj()[0]
    vmilieu = traj.getTraj()[0] - traj.getTraj()[-1]
    vfin = traj.getTraj()[-1] - traj.getTraj()[-2]

    # On calcule la similarité cosinus entre les vecteurs
    dot1 = (vdebut.dot(vmilieu)) / (vdebut.magnitude * vmilieu.magnitude)
    dot2 = (vmilieu.dot(vfin)) / (vmilieu.magnitude * vfin.magnitude)

    print("\n-----------------------------------------------------------------------------------------")
    print("La distance entre le centre du dinucléotide du début et celui de fin du meilleur individu : ",
          (traj.getTraj()[0] - traj.getTraj()[-1]).magnitude)
    print(f"La similarité cosinus :\nEntre v_fin et v_milieu : {dot1}\nEntre v_milieu et v_debut : {dot2}")
    print("-----------------------------------------------------------------------------------------")

    # On affiche la trajectoire
    traj.draw()
    traj.write(args.filename+".png")


if __name__ == "__main__":
    main()
