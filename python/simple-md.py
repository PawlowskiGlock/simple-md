#!/usr/bin/env python3
"""
@author: micabo
@date: 28.03.17
"""

import numpy as np

# atom / fragment mass
m = 15.035

# bond constants
k = 3e5
b0 = 0.153

# lennard-jones constants
C12 = 33.570436e-6
C6 = 0.0099161764


class Atom(object):
    id = 0

    def __init__(self, mass, pos, vel):
        Atom.id += 1
        self.id = Atom.id
        self.mass = float(mass)
        self.pos = np.array(pos, dtype=float)
        self.vel = np.array(vel, dtype=float)
        self.force = np.array((0.0, 0.0, 0.0), dtype=float)

    def __del__(self):
        Atom.id -= 1


class Bond(object):

    def __init__(self, k, b0, atom1, atom2):
        self.k = k
        self.b0 = b0
        self.atom1 = atom1
        self.atom2 = atom2


class Molecule(object):

    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds


class System(object):

    def __init__(self, mols):
        self.mols = mols

    def _bond_force(self, bond):
        r = bond.atom1.pos - bond.atom2.pos
        b = np.sqrt(r.dot(r))
        r /= b  # normalize
        F = bond.k * (b - bond.b0) * r
        bond.atom1.force -= F
        bond.atom2.force += F

    def _lj_force(self, atom1, atom2):
        dist = atom1.pos - atom2.pos
        r = np.sqrt(dist.dot(dist))
        dist /= r  # normalize
        C = (12*C12/r**14 - 6*C6/r**8) * dist
        atom1.force += C
        atom2.force -= C

    def _kick(self, dt):
        for i in range(len(self.mols)):
            # bonded interaction
            for bond in self.mols[i].bonds:
                self._bond_force(bond)  # includes side-effect
            # LJ-interaction - only INTERmolecular
            for j in range(i+1, len(self.mols)):
                for atom1 in self.mols[i].atoms:
                    for atom2 in self.mols[j].atoms:
                        self._lj_force(atom1, atom2)  # includes side-effect
            # propagate force into velocity
            for atom in self.mols[i].atoms:
                atom.vel += atom.force * dt / atom.mass
                atom.force.fill(0.0)  # reset force to 0

    def _drift(self, dt):
        for mol in self.mols:
            for atom in mol.atoms:
                atom.pos += atom.vel * dt

    def _bond_energy_x2(self, bond):
        r = bond.atom1.pos - bond.atom2.pos
        b = np.sqrt(r.dot(r))
        E = bond.k * (b - bond.b0)**2
        return E

    def _lj_energy(self, atom1, atom2):
        r = atom1.pos - atom2.pos
        r = np.sqrt(r.dot(r))
        E = (12*C12/r**12 - C6/r**6)
        return E

    def kinetic_energy(self):
        k_energ = 0.0
        for mol in self.mols:
            for atom in mol.atoms:
                k_energ += atom.mass * atom.vel.dot(atom.vel)
        return k_energ / 2.0

    def potential_energy(self):
        pot_energ = 0.0
        for i in range(len(self.mols)):
            # bonded energy
            for bond in self.mols[i].bonds:
                pot_energ += self._bond_energy_x2(bond)
            pot_energ /= 2.0
            # LJ-energy
            for j in range(i+1, len(self.mols)):
                for atom1 in self.mols[i].atoms:
                    for atom2 in self.mols[j].atoms:
                        pot_energ += self._lj_energy(atom1, atom2)
        return pot_energ

    def print_atoms(self, step, timestep):
        print(Atom.id)
        # think about units!
        energy = self.kinetic_energy() + self.potential_energy()
        ite = (step, step*timestep, energy)  # i, time, energy
        print("i = {:10d}, time = {:10.4f}, energy = {:10.6f}".format(*ite))
        for mol in self.mols:
            for atom in mol.atoms:
                pos = 10 * atom.pos  # nm -> angstrom
                print("CH3 {:+14.10f} {:+14.10f} {:+14.10f}".format(*pos))

    def md(self, nsteps, timestep):
        timestep_2 = timestep / 2.0
        for step in range(1, nsteps+1):
            self._drift(timestep_2)
            self._kick(timestep)
            self._drift(timestep_2)
            self.print_atoms(step, timestep)


if __name__ == "__main__":
    # create atoms
    """
    # scenario 1
    a1 = Atom(m, ( 0.0000,  0.0000,  0.0000), ( 0.0000,  0.0000,  0.0000))
    a2 = Atom(m, ( 0.1530,  0.0000,  0.0000), ( 0.0000,  0.0000,  0.0000))
    a3 = Atom(m, ( 2.0000,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))
    a4 = Atom(m, ( 2.1530,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))

    # scenario 2
    a1 = Atom(m, ( 0.0000, -0.0765,  0.0000), ( 0.0000,  0.0000,  0.0000))
    a2 = Atom(m, ( 0.0000,  0.0765,  0.0000), ( 0.0000,  0.0000,  0.0000))
    a3 = Atom(m, ( 2.0000,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))
    a4 = Atom(m, ( 2.1530,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))

    #scenario 3
    a1 = Atom(m, ( 0.0000,  0.0000,  0.0000), ( 0.0000,  0.0000,  0.0000))
    a2 = Atom(m, ( 0.0000,  0.1530,  0.0000), ( 0.0000,  0.0000,  0.0000))
    a3 = Atom(m, ( 2.0000,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))
    a4 = Atom(m, ( 2.1530,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))
    """
    #scenario 4
    a1 = Atom(m, ( 0.0000,  0.0000,  0.0000), ( 0.0000,  0.5000,  0.0000))
    a2 = Atom(m, ( 0.0000,  0.1530,  0.0000), ( 0.0000, -0.5000,  0.0000))
    a3 = Atom(m, ( 2.0000,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))
    a4 = Atom(m, ( 2.1530,  0.0000,  0.0000), (-1.0000,  0.0000,  0.0000))

    # create bonds
    b12 = Bond(k, b0, a1, a2)
    b34 = Bond(k, b0, a3, a4)

    # create molecules
    m1 = Molecule([a1, a2], [b12])
    m2 = Molecule([a3, a4], [b34])

    # create system
    sim = System([m1, m2])
    sim.md(2000, 0.002)
