import os
import sys
from math import sqrt, log
import matplotlib.pyplot as plt

CUTOFF = 15
EPS = 1e-5

def read_pdb(fname):
    """
    Return a list of (resname, atomname, resnum, (x, y, z))
       for all 'ATOM' entries in FNAME.

    Note that in a research project, I would generally use a library
    like BioPython or rdkit to read pdb files, but this simple
    parser is sufficient for our purpose here and avoids having to
    install additional libraries.
    """
    coords = []
    with open(fname) as fp:
        for line in fp:
            # PDB files have fixed width columns!
            if len(line) < 4 or line[:4] != 'ATOM': continue
            atomname = line[12:16].strip()
            resname = line[17:20].strip()
            resnum  = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords += [(resname, atomname, resnum, (x, y, z))]
    return coords

def read_train_set():
    """
    Parse all structures in the training set. Returns a list
    of the output from read_pdb for each structure.
    """
    root = 'train-pdbs'
    assert os.path.exists(root)
    structures = []
    for pdb_fname in os.listdir(root):
        if pdb_fname[0] == '.': continue
        structures += [read_pdb(f'{root}/{pdb_fname}')]
    return structures

def read_test_set():
    """
    Parse native and decoy structures in the test set.
    Returns a dict with entry "native" giving the native
    structure and entry "decoys" giving a list of decoys.
    """
    root = 'test-pdbs-moulder'
    test_set = {}
    for pdb in os.listdir(root):
        if pdb[0] == 'l': continue
        test_set[pdb] = {
            'structs': [read_pdb(f'{root}/{pdb}/model{i}.pdb')
                        for i in range(2, 300)]
        }
    # Read RMSDs.
    with open(f'{root}/list-without-native') as fp:
        for line in fp:
            if line[0] == '*':
                pdb = line[1:].strip()
                test_set[pdb]['rmsds'] = []
            else:
                # One RMSD is missing, for some reason...
                rmsd = line.split('|')[1].strip()
                rmsd = float(rmsd) if rmsd else 10.0
                if len(test_set[pdb]['rmsds']) < len(test_set[pdb]['structs']):
                    test_set[pdb]['rmsds'] += [rmsd]
    return test_set

def get_coords(resname, structure):
    """
    Return a list of (x, y, z) coordinates occurances of the
    given resname in the given structure.
    """
    # Glycine doesn't have a CB.
    if resname == 'GLY':
        atomname = 'CA'
    else:
        atomname = 'CB'

    coords = []
    for _resname, _atomname, resnum, coord in structure:
        if _resname == resname and  _atomname == atomname:
            coords += [(resnum, coord)]
    return coords

def get_distances(coords1, coords2):
    """
    Return a list of all pairwise distances between coordinates in coords1
    and coords2. Exclude neighoring residues.
    """
    dists = []
    for resnum1, coord1 in coords1:
       for resnum2, coord2 in coords2:
           if resnum1 + 1 < resnum2 or resnum2 < resnum1 - 1:
               dists += [sum((coord1[i]-coord2[i])**2 for i in range(3))**0.5]
    return dists

def histogram(distances):
    """
    Returns a histogram of values where each range spans 1 unit.
    Values greater than cutoff are ignored.
    Entry y[i] corresponds to the number of values between x[i] and x[i+1]
    """
    x = list(range(CUTOFF))
    y = [EPS]*len(x) # Start with a small count to avoid divide by zeros errors.
    for distance in distances:
        if int(distance) < CUTOFF:
            y[int(distance)] += 1
    return x, y

#################################################################

def plot_pair(title, x, y_obs, y_exp, energy):
    """
    Plot observed and expected in top panel and energy in
    bottom panel.
    """
    if not y_exp:
        print('Expected distribution not implemented')
        y_exp = [0]*len(y_obs)

    if not energy:
        print('Energy not implemented')
        energy = [0]*len(y_obs)

    f, ax = plt.subplots(2, sharex=True)
    ax[0].plot(x, y_obs, label='observed')
    ax[0].plot(x, y_exp, label='expected')
    ax[0].set_ylabel('Count')
    ax[0].legend(loc='best')

    ax[1].plot(x, energy)
    ax[1].plot(x, [0 for _ in energy], ls='--', c='k')
    ax[1].set_ylabel('Energy (Arbitrary Units)')
    ax[1].set_xlabel('Distance (Angstroms)')
    plt.ylim(-3, 3)
    ax[0].set_title(title)
    plt.show()

def plot_all(energies):
    """
    Plot energy values for many resname pairs.
    energies should be a dictionary with keys (resname1, resname2)
    and values indicating the energy as a function of distance.
    """
    resnames = sorted(set(r for r in energies for r in r))

    f, ax = plt.subplots(len(resnames), len(resnames), figsize=(10, 10), sharex=True, sharey=True)

    for i, resname1 in enumerate(resnames):
        for j, resname2 in enumerate(resnames):
            if i >= j:
                energy = energies[(resname1, resname2)]
                ax[i, j].plot(energy)
                ax[i, j].plot([0 for i in energy], ls='--', c='k')

    for i, resname in enumerate(resnames):
        ax[i, 0].set_ylabel(resname)
        ax[-1, i].set_xlabel(resname)

    plt.ylim(-2, 2)
    plt.show()

def write_energies(energies):
    with open('energies.txt', 'w') as fp:
        for (resname1, resname2), energy in energies.items():
            energy = ','.join(map(str, energy))
            fp.write(f'{resname1},{resname2}:{energy}\n')

def read_energies():
    energies = {}
    with open('energies.txt') as fp:
        for line in fp:
            resname, energy = line.strip().split(':')
            resname1, resname2 = resname.split(',')
            energy = [float(x) for x in energy.split(',')]
            energies[(resname1, resname2)] = energy
    return energies

def score(energies):
    """
    For each entry in the test set, compute the rank of the native
    structure amongst the decoys, i.e. 1 + the number of decoy
    structures with a more favorable score than the native structure.
    """
    test_cases = read_test_set()
    f, ax = plt.subplots(4, 5, figsize=(12, 8))
    for i, pdb in enumerate(test_cases):
        scores = []
        best_score = float('inf')
        best_idx = None
        for j, structure in enumerate(test_cases[pdb]['structs']):
            score = scorer(structure, energies)
            scores += [score]
            if score < best_score:
                best_score = score
                best_idx = j

        pred_rmsd = test_cases[pdb]['rmsds'][best_idx]
        min_rmsd = min(test_cases[pdb]['rmsds'])
        median_rmsd = sorted(test_cases[pdb]['rmsds'])[int(len(test_cases[pdb]['rmsds'])/2)]
        print(f'PDB {pdb}: prediction rmsd = {pred_rmsd}, best possible rmsd = {min_rmsd}, median rmsd = {median_rmsd}')

        _ax = ax[i%4, int(i/4)]
        _ax.scatter(test_cases[pdb]['rmsds'], scores)
        _ax.set_title(pdb)
    for i in range(0, 5):
        ax[-1, i].set_xlabel('rmsd')
    for i in range(0, 4):
        ax[i, 0].set_ylabel('energy')
    plt.show()

#####################################################################

def observed_distribution(resname1, resname2, train_structures):
    """
    Compute a histogram of the observed distances between each pair
    of residues with residue type resname1 and resname2 for each
    structure in the training set.
    """

    ##################################################################
    # The get_coords, get_distances, and histogram function will be useful.

    distances = []
    # Edit here
    for structure in train_structures:
        coords_resname1 = get_coords(resname1, structure)
        coords_resname2 = get_coords(resname2, structure)
        dist_list = get_distances(coords_resname1, coords_resname2)
        # print("Dist List is :")
        # print(dist_list)
        for dist in dist_list:
            # print(dist)
            distances.append(dist)
    x, y_obs = histogram(distances)
    ##################################################################
    return x, y_obs # The output from `histogram`

def expected_distribution(x):
    """
    Given x a list of integers, return a list of the  expected number
    of occurances between x[i] and x[i] + 1.

    Ignore the 4 * pi multiplier.
    """
    #############################################################
    y_exp = []
    for i in x:
        # Edit here
        pass # this is a placeholder; remove this line once you start implementing your function
    #############################################################
    return y_exp

def distributions_to_energy(y_obs, y_exp):
    """
    Given the observed and expected counts, return a list of the energy
    for each distance bin.
    """
    #############################################################
    # The log function will be useful.
    energy = []
    # Edit here

    #############################################################
    return energy

def pair_energy(resname1, resname2):
    """
    Compute the observed and expected distributions and corresponding
    energies for resname1 and resname2.
    """
    train_structures = read_train_set()
    x, y_obs = observed_distribution(resname1, resname2, train_structures)
    y_exp = expected_distribution(x)
    y_exp = [y*y_obs[-1]/y_exp[-1] for y in y_exp]
    energy = distributions_to_energy(y_obs, y_exp)
    return x, y_obs, y_exp, energy

def all_pair_energies():
    """
    Compute the energies for all pairs of residue types.
    Returns a dict with keys (resname1, resname2) mapping
    to a list of energies for each distance bin.
    """
    resnames = ['LYS', 'ARG', 'ASP', 'GLU', 'TRP',
                'PHE', 'TYR', 'LEU', 'ILE', 'MET',
                'VAL', 'PRO', 'GLY', 'CYS', 'ALA',
                'SER', 'ASN', 'GLN', 'THR', 'HIS']
    energies = {}
    for resname1 in resnames:
        for resname2 in resnames:
            print(resname1, resname2)
            if resname2 < resname1:
                continue
            x, y_obs, y_exp, energy = pair_energy(resname1, resname2)
            energies[(resname1, resname2)] = energy
            energies[(resname2, resname1)] = energy
    return energies

def scorer(structure, energies):
    """
    Given a structure and energies (as returned by all_pair_energies), return
    the score for the structure

       E = \sum_i \sum_j pair_energy(distance_ij, resname_i, resname_j)

    Note that it will be easier to write your code using the equivalent formulation

       E = \sum_{resname1} \sum_{resname2}
               \sum_{i where resname = resname1} \sum_{j where resname = resname2}
                      pair_energy(distance_ij, resname_i, resname_j)
    """

    ###############################################################
    # The get_coords and get_distances function will be helpful.
    # Part of your implementation should be similar to your implementation of
    # observed_distribution.
    # To lookup an energy for a given distance, you can convert the
    # distance to an index using i = int(distance).
    # Distances greater than or equal to the cutoff should be ignored!
    score = 0
    for (resname1, resname2), energy in energies.items():
        if resname1 <= resname2:
            # Edit here
            pass # remove this line once you start implementing your function

    ###############################################################
    return score

if __name__ == '__main__':
    if len(sys.argv) == 3:
        resname1 = sys.argv[1]
        resname2 = sys.argv[2]
        x, y_obs, y_exp, energy = pair_energy(resname1, resname2)
        plot_pair(f'{resname1}-{resname2}', x, y_obs, y_exp, energy)
    elif len(sys.argv) == 2 and sys.argv[1] == 'compute':
        energies = all_pair_energies()
        write_energies(energies)
        plot_all(energies)
    elif len(sys.argv) == 2 and sys.argv[1] == 'score':
        energies = read_energies()
        score(energies)
    else:
        print('Invalid arguments.')
