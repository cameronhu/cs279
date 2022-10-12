from pymol import cmd
import matplotlib.pyplot as plt

def get_residue_numbers(sel):
    """
    Return a list of integers giving residue numbers for all protein
    residues in the given selection.
    """
    resnums = {'resnums': []}
    cmd.iterate(f'{sel} and name CA', 'resnums.append(resi)', space=resnums)
    resnums = resnums['resnums']
    return [int(resnum) for resnum in resnums]

def get_n_states(sel):
    """
    Returns the number of states for the given selection.
    """
    return cmd.count_states(sel)

def plot_distance(entry, sel1, sel2):
    n_states = get_n_states(entry)
    cmd.distance(f'{entry} and {sel1}', f'{entry} and {sel2}')
    distances = []
    for state in range(1, n_states+1):
        distances += [cmd.get_distance(f'{entry} and {sel1}', f'{entry} and {sel2}', state=state)]

    plt.plot(distances)
    plt.ylabel(f'Distance from {sel1} to {sel2} (Angstroms)')
    plt.xlabel('Time')
    plt.title(entry)
    plt.show()

def plot_rmsd(ref_sel, query_sel):
    n_states = get_n_states(query_sel)
    rmsds = []
    for	state in range(1, n_states+1):
        rmsds += [cmd.align(query_sel, ref_sel, target_state=0, mobile_state=state, cycles=0, transform=0)[0]]

    plt.plot(rmsds)
    plt.ylabel(f'RMSD (Angstroms)')
    plt.xlabel('Time')
    plt.show()

def ramachandran(sel):
    """
    Produce a Ramachandran plot for residues in the given selection.
    """
    cmd.delete('phi')
    cmd.delete('psi')
    resnums = get_residue_numbers(sel)
    phis, psis = [], []

    # For each residue in resnums, add the phi and psi angle to phis and psis, respectively.

    # PyMol has two commands related to dihedral angles:
    # cmd.dihedral(name, sel1, sel2, sel3, sel4) will plot the dihedral on the
    #     structure however it will not return the value.
    # cmd.get_dihedral(sel1, sel2, sel3, sel4) will return the dihedral but will
    #     not show it on the structure.

    # Only cmd.get_dihedral is strictly required in your implementation, but
    # we highly recommend that you call both commands with the same selections
    # so that you can visually see the angles for debugging purposes.
    # Note that cmd.dihedral has an additional ``name'' argument, which you
    # should set to "phi" for the phi angles and psi for the "psi" angles.

    # Some tips on what various error messages mean:
    # "Error: Selection 1: Not found": The first selection matches no atoms.
    # "Error: Selection 1: Invalid selection name": The first selection matches multiple atoms.
    # Equivalent messages for Selection 2 mean the second selection is invalid, and so on.

    ############################################################################
    # Edit here.

    #psi angles measure the dihedral angle between c(i-1), n(i), calpha(i), c(i)
    #psi angles measure the dihedral angle between n(i), calpha(i), c(i), n(i+1)

    for resnum in resnums:
        #i-1 acyl carbon
        cmd.select(f'{resnum - 1}c', f'residue {resnum - 1} and name c')
        #nitrogen
        cmd.select(f'{resnum}n', f'residue {resnum} and name n')
        #calpha
        cmd.select(f'{resnum}ca', f'residue {resnum} and name ca')
        #acyl_carbon
        cmd.select(f'{resnum}c', f'residue {resnum} and name c')
        #i+1 nitrogen
        cmd.select(f'{resnum + 1}n', f'residue {resnum + 1} and name n')
        phi_angle = cmd.get_dihedral(f'{resnum-1}c', f'{resnum}n', f'{resnum}ca', f'{resnum}c')
        cmd.dihedral(f'{resnum}_phi', f'{resnum-1}c', f'{resnum}n', f'{resnum}ca', f'{resnum}c')
        psi_angle = cmd.get_dihedral(f'{resnum}n', f'{resnum}ca', f'{resnum}c', f'{resnum+1}n')
        cmd.dihedral(f'{resnum}_psi', f'{resnum}n', f'{resnum}ca', f'{resnum}c', f'{resnum+1}n')
        phis.append(phi_angle)
        psis.append(psi_angle)

    ############################################################################
    plt.scatter(phis, psis)
    plt.xlabel('phi')
    plt.ylabel('psi')
    plt.ylim(-180, 180)
    plt.xlim(-180, 180)
    plt.gca().set_aspect('equal')
    plt.show()

cmd.extend('plot_distance', plot_distance)
cmd.extend('plot_rmsd', plot_rmsd)
cmd.extend('ramachandran', ramachandran)
