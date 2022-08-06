from ase import io
atoms = io.read('input.cif')
atoms.write('INPUT_STRUCTURE', format = 'vasp')
