from pymol import cmd
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from numpy import random

def compute_sasa(pdbs_folder):
    pdb_files = [f for f in listdir(pdbs_folder) if isfile(join(pdbs_folder, f))]
    amino_acids = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN',
    'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    aa_mw = [174.2, 155.2, 146.2, 133.1, 147.1, 105.1, 119.1, 132.1, 146.2, 121.2, 75.1,
    115.1, 89.1, 117.1, 131.2, 131.2, 149.2, 165.2, 181.2, 204.2]
    aa_MaxASA = []
    polar_aa = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU']
    sasa_dict = {}
    for amino_acid in amino_acids:
        sasa_dict[amino_acid] = [0, 0, 0]

    for file in pdb_files:
        cmd.load(pdbs_folder +'/' + file)
        #need to set dot_solvent to 1 for SASA, 0 is total available molecular surface area
        cmd.set('dot_solvent', '1')
        #cmd.set('dot_density', '1')
        cmd.h_add()
        for aa in amino_acids:
            cmd.select(f'{aa}_sel', f'resn {aa}')
            num_residue = cmd.count_atoms(f'{aa}_sel and name ca')
            sasa_dict[aa][0] += num_residue
            sasa_dict[aa][1] += cmd.get_area(f'{aa}_sel')
        cmd.delete('all')

    sasa_averages = {}
    colors = []
    for amino_acid in amino_acids:
        sasa_averages[amino_acid] = [] #For each amino acid, store non_normalized SASA, normalized SASA by MW, and normalized SASA by Max ASA
        this_sasa = sasa_dict[amino_acid][1] / sasa_dict[amino_acid][0]
        sasa_averages[amino_acid].append(this_sasa)
        if amino_acid in polar_aa:
            colors.append('orange')
        else:
            colors.append('blue')

    for key in sasa_averages.get_keys():
        sasa_averages_norm_mw.append(sasa_averages[i] / aa_mw[i])
        sasa_averages_norm_total_sasa.append(i)

    #Random values to lists for matplotlib testing
    # sasa_averages = []
    # sasa_average_norm = []
    # colors = []
    # for i in range(20):
    #     sasa_averages.append(random.randint(100))
    #     sasa_average_norm.append(random.randint(100))
    #     colors.append('purple')

    fig1, ax1 = plt.subplots(figsize = (16, 9))
    fig2, ax2 = plt.subplots(figsize = (16, 9))
    fig3, ax3 = plt.subplots(figsize = (16, 9))
    ax1.set_title("Average SASA of Different Amino Acid Residues Across the Given PDBs")
    ax1.bar(amino_acids, sasa_averages, color=colors)
    ax2.set_title("Average SASA Normalized by Molecular Weight of Different Amino Acid Residues Across the Given PDBs")
    ax2.bar(amino_acids, sasa_averages_norm_mw, color = colors)
    ax3.set_title("Average SASA Normalized by Total Accessible Surface Area of Different Amino Acid Residues Across the Given PDBs")
    ax3.bar(amino_acids, sasa_averages_norm_total_sasa, color = colors)
    set_label_titles(ax1)
    set_label_titles(ax2)
    set_label_titles(ax3)

    plt.show()

def set_label_titles(axis):
    axis.xaxis.set_tick_params(pad = 5)
    axis.set_xlabel("Amino Acid Residue")
    axis.set_ylabel("Average SASA (Square Angstroms)")

cmd.extend('compute_sasa', compute_sasa)
