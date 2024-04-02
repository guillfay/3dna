# ST2 (Théorie des Jeux) - EI Algorithmique Génétique

![CentraleSupelec Logo](https://www.centralesupelec.fr/sites/all/themes/cs_theme/medias/common/images/intro/logo_nouveau.jpg)

- Le rapport du projet est disponible :point_right: [ici](Rapport.pdf)
- Le support de la soutenance est disponible :pôint_right: [ici](Diapo_soutenance.pdf)

## Initialisation

- Créer un environnement virtuel Python
```sh
python -m venv venv
```
- Activer l'environnement
    - Windows
    ```powershell
    .\venv\Scripts\activate
    ```
    - Linux / macos
    ```sh
    source ./venv/Scripts/activate
    ```

- Installer les librairies Python
```sh
pip install -r requirements.txt
```

## Recherche d'une solution

- Algorithme du recuit simulé
```sh
python 3dna <répertoire du fichier .fasta> recuit
```
Exemple
```sh
python 3dna ./data/plasmid_8k.fasta recuit
```

- Algorithme génétique
```sh
python 3dna <répertoire du fichier .fasta> genetique <stratégie de sélection>
```
Liste des stratégies de sélection
```
elitisme (défaut)
roulette
rang
tournoi
```

Exemple
```sh
python 3dna ./data/plasmid_8k.fasta genetique tournoi
```

## Valeur retournée

A la fin de l'exécution d'un algorithme, plusieurs valeurs sont imprimées :
- La table des rotations 3D en format `json`
- La différence entre les valeurs de *Twist* et *Wedge* et leur seuil maximal pour chaque paire de nucléotides (La valeur de *Direction* est constante)
- La distance entre le centre du dinucléotide du début et celui de fin.
- En notant
    - $p_i$ la position du centre du ième dinucléotide
    - $n$ le nombre de nucléotides dans la séquence
    - $v_{fin} = p_{n - 1} - p_{n - 2}$
    - $v_{milieu} = p_0 - p_{n - 1}$
    - $v_{debut} = p_1 - p_0$

    La similarité cosinus entre les vecteurs $v_{fin}$ et $v_{milieu}$, et $v_{milieu}$ et $v_{debut}$

## Configuration des algorithmes
Les valeurs internes utilisées pour l'algorithme peuvent être modifiés dans le fichier `3dna/__main__.py`.
- Algorithme de recuit simulé
```python
recuit = Recuit(
    seq=seq,
    k_max=1000,
    e_min=0.1,
    temp_init=300,
    refroidissement=0.99,
    dist_min=0.,
    relier=1
)
```

|argument|type|description|
|-|-|-|
|seq|str|La séquence de nucléotides|
|k_max|int|Le nombre d'itération maximale|
|e_min|float|Le seuil d'énergie|
|temp_init|float|La température initiale|
|refroidissement|float|Le coefficient de refroidissement|
|dist_min|float|Le seuil de distance|
|Relier|int|Le nombre de nucléotides reliés entre la fin et le début de la séquence|

- Algorithme génétique
```python
genetique = Genetique(
    seq=seq,
    methode_utilisee=args.methode_utilisee,
    nbr_generation_max=150,
    N=50,
    probabilite_mutation_initiale=0.05,
    probabilite_mutation_finale=0.05,
    relier = 1
)
```

|argument|type|description|
|-|-|-|
|seq|str|La séquence de nucléotides|
|methode_utilisee|str|Le nom de la stratégie de sélection|
|nbr_generation_max|int|Le nombre de génération maximal|
|N|int|La taille de la population|
|probabilite_mutation_initiale|float|La probabilité de mutation initiale|
|probabilite_mutation_finale|float|La probabilité de mutation finale|
|Relier|int|Le nombre de nucléotides reliés entre la fin et le début de la séquence|

*Note : La probabilité de mutation évolue de manière linéaire entre la valeur initiale et la valeur finale*

## Test

Un test a été implémenté grâce à la librairie `pytest`.

Il peut être exécuté avec

```
cd 3dna
pytest
```
Un rapport de couverture sous format `html` peut être généré avec
```
python -m pytest --cov=. --cov-report html
```

## Statistiques

Les algorithmes peuvent être testés avec plusieurs centaines d'itérations afin d'obtenir une distribution de la distance entre le centre du premier et du dernier dinucléotide, et du minimum de la similarité cosinus.

```bash
python 3dna/stat.py
python 3dna/plot.py
```
