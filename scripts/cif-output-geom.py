from ase import io
atoms = io.read('POSCAR')
atoms.write('OUTPUT_STRUCTURE', format = 'cif')
