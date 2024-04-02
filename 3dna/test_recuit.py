from RotTable import RotTable
from Traj3D import Traj3D
from recuit import Recuit
from random import seed
from math import exp

# Seed pour la reproductibilité des résultats de test
seed(91)

def test_energie():
    # Initialisation des paramètres
    seq=''.join([line.rstrip('\n') for line in open(r'..\data\plasmid_8k.fasta')][1:])
    k_max=1
    e_min=1
    temp_init=1
    refroidissement=1
    dist_min=1
    
    # Test n°1
    test_recuit=Recuit(seq, k_max, e_min, temp_init, refroidissement,dist_min)
    assert test_recuit.energie(RotTable())==(5987.480078776554, 5990.480078776554)

    dist_min=0
    # Test n°2
    test_recuit=Recuit(seq, k_max, e_min, temp_init, refroidissement,dist_min)
    assert test_recuit.energie(RotTable())==(5987.480078776554, 5990.480078776554)
   
def test_voisin():
    # Initialisation des paramètres
    seq=''.join([line.rstrip('\n') for line in open(r'..\data\plasmid_8k.fasta')][1:])
    k_max=1
    e_min=1
    temp_init=1
    refroidissement=1
    dist_min=1
    relier=1

    # Test
    test_recuit=Recuit(seq, k_max, e_min, temp_init, refroidissement,dist_min,relier)
    assert test_recuit.voisin(RotTable(),100).rot_table=={'AA': [35.68, 6.6000000000000005, -154, 0.06, 0.6, 0], 'AC': [35.699999999999996, -3.9, 143, 1.3, 5, 0], 'AG': [29.2, 5.4, 2, 1.5, 3, 0], 'AT': [32.6, 4.6, 0, 1.1, 2, 0], 'CA': [35.4, 37.5, -64, 0.9, 34, 0], 'CC': [33.6, 0.0, -57, 0.07, 2.1, 0], 'CG': [28.7, 5.2, 0, 1.1, 1.5, 0], 'CT': [29.2, 5.4, -2, 1.5, 3, 0], 'GA': [37.8, 11.3, 120, 0.9, 6, 0], 'GC': [41.2, 6.275, 180, 1.2, 1.275, 0], 'GG': [33.6, 0.0, 57, 0.07, 2.1, 0], 'GT': [35.699999999999996, -3.9, -143, 1.3, 5, 0], 'TA': [37.1, -1.1, 0, 1.1, 2, 0], 'TC': [37.8, 11.3, -120, 0.9, 6, 0], 'TG': [35.4, 37.5, 64, 0.9, 34, 0], 'TT': [35.68, 6.6000000000000005, 154, 0.06, 0.6, 0]}

def test_P():
    # Initialisation des paramètres
    seq=''.join([line.rstrip('\n') for line in open(r'..\data\plasmid_8k.fasta')][1:])
    k_max=1
    e_min=1
    temp_init=1
    refroidissement=1
    dist_min=1
    relier=1

    # Test
    test_recuit=Recuit(seq, k_max, e_min, temp_init, refroidissement,dist_min,relier)
    assert test_recuit.P(-1,2)==1
    assert test_recuit.P(1,2)==exp(-1/2)

def test_executer():
    # Initialisation des paramètres
    seq=''.join([line.rstrip('\n') for line in open(r'..\data\plasmid_8k.fasta')][1:])
    k_max=1
    e_min=1
    temp_init=1
    refroidissement=1
    dist_min=1
    relier=1

    # Test n°1
    test_recuit=Recuit(seq, k_max, e_min, temp_init, refroidissement,dist_min,relier)
    assert test_recuit.executer().rot_table=={'AA': [35.62, 7.2, -154, 0.06, 0.6, 0], 'AC': [34.4, 1.1, 143, 1.3, 5, 0], 'AG': [27.7, 8.4, 2, 1.5, 3, 0], 'AT': [31.5, 2.6, 0, 1.1, 2, 0], 'CA': [34.5, 3.5, -64, 0.9, 34, 0], 'CC': [33.67, 2.1, -57, 0.07, 2.1, 0], 'CG': [29.8, 6.7, 0, 1.1, 1.5, 0], 'CT': [27.7, 8.4, -2, 1.5, 3, 0], 'GA': [36.9, 5.3, 120, 0.9, 6, 0], 'GC': [40, 5, 180, 1.2, 1.275, 0], 'GG': [33.67, 2.1, 57, 0.07, 2.1, 0], 'GT': [34.4, 1.1, -143, 1.3, 5, 0], 'TA': [36, 0.9, 0, 1.1, 2, 0], 'TC': [36.9, 5.3, -120, 0.9, 6, 0], 'TG': [34.5, 3.5, 64, 0.9, 34, 0], 'TT': [35.62, 7.2, 154, 0.06, 0.6, 0]}
    
    temp_init=10000
    # Test n°2
    test_recuit=Recuit(seq, k_max, e_min, temp_init, refroidissement,dist_min,relier)
    assert test_recuit.executer().rot_table=={'AA': [35.6258723685025, 7.133218258737684, -154, 0.06, 0.6, 0], 'AC': [34.856033195455424, 1.5526980994848585, 143, 1.3, 5, 0], 'AG': [28.415129954458383, 8.863504956186638, 2, 1.5, 3, 0], 'AT': [31.473929741161417, 2.160910617661816, 0, 1.1, 2, 0], 'CA': [34.82516341165077, 17.632652071392595, -64, 0.9, 34, 0], 'CC': [33.64726924542647, 2.3270116148355133, -57, 0.07, 2.1, 0], 'CG': [29.85747152236308, 6.601506677853863, 0, 1.1, 1.5, 0], 'CT': [28.415129954458383, 8.863504956186638, -2, 1.5, 3, 0], 'GA': [37.18593497127491, 4.123281251727021, 120, 0.9, 6, 0], 'GC': [40.25935916389681, 4.635416200637767, 180, 1.2, 1.275, 0], 'GG': [33.64726924542647, 2.3270116148355133, 57, 0.07, 2.1, 0], 'GT': [34.856033195455424, 1.5526980994848585, -143, 1.3, 5, 0], 'TA': [35.15849209286032, 2.025997573459315, 0, 1.1, 2, 0], 'TC': [37.18593497127491, 4.123281251727021, -120, 0.9, 6, 0], 'TG': [34.82516341165077, 17.632652071392595, 64, 0.9, 34, 0], 'TT': [35.6258723685025, 7.133218258737684, 154, 0.06, 0.6, 0]}
    # {'AA': [35.62530865041074, 7.319907819695622, -154, 0.06, 0.6, 0], 'AC': [34.82812244838177, 0.4630498229477823, 143, 1.3, 5, 0], 'AG': [27.07788607359709, 7.011380234733055, 2, 1.5, 3, 0], 'AT': [30.835663425339746, 2.4201677565617152, 0, 1.1, 2, 0], 'CA': [34.55279996109801, -11.659571432015275, -64, 0.9, 34, 0], 'CC': [33.625542978007545, 1.537683929641682, -57, 0.07, 2.1, 0], 'CG': [29.650451940001584, 6.295382019308093, 0, 1.1, 1.5, 0], 'CT': [27.07788607359709, 7.011380234733055, -2, 1.5, 3, 0], 'GA': [37.3449537691803, 3.132160446333805, 120, 0.9, 6, 0], 'GC': [40.5597570978094, 4.812761178018319, 180, 1.2, 1.275, 0], 'GG': [33.625542978007545, 1.537683929641682, 57, 0.07, 2.1, 0], 'GT': [34.82812244838177, 0.4630498229477823, -143, 1.3, 5, 0], 'TA': [35.68142599176165, 1.110505173032761, 0, 1.1, 2, 0], 'TC': [37.3449537691803, 3.132160446333805, -120, 0.9, 6, 0], 'TG': [34.55279996109801, -11.659571432015275, 64, 0.9, 34, 0], 'TT': [35.62530865041074, 7.319907819695622, 154, 0.06, 0.6, 0]}