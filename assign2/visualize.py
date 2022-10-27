from pymol import cmd
import sys
sys.path.append('.')
from statistical_potential import read_energies

energies = read_energies()

def get_atoms(sel):
    """
    Return a list of "entry and chain/resi/name" for all atoms
    in the given selection.
    """
    atoms = {'resnums': []}
    cmd.iterate(f'{sel} and (name CB or (resname GLY and name CA))',
                'resnums.append((resn, resi+"/"+name))', space=atoms)
    atoms = atoms['resnums']
    return atoms

def show_energies(entry):
    scores = []
    atoms = get_atoms(entry)
    for i, (resname1, atom1) in enumerate(atoms):
        for resname2, atom2 in atoms[i+2:]:
            distance = cmd.get_distance(f'{entry} and {atom1}', f'{entry} and {atom2}')
            print(distance, resname1, resname2)
            if (resname1, resname2) in energies and distance < 15:
                energy = energies[(resname1, resname2)][int(distance)]
                print(energy)

                if energy < -0.75:
                    cmd.distance(f'{entry}_favorable', f'{entry} and {atom1}', f'{entry} and {atom2}')
                    cmd.show('sticks', f'br. {entry} and {atom1}')
                    cmd.show('sticks', f'br. {entry} and {atom2}')

    cmd.color('green', f'{entry}_favorable')

cmd.extend('show_energies', show_energies)
cmd.auto_arg[0]['show_energies'] = [ cmd.object_sc, 'object', ''] # autocomplete
