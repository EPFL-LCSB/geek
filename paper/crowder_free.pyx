from __future__ import print_function

import numpy as np
import random

import numpy.random as rnd
from math import sqrt, exp, log, floor
import copy

AVOGADRO_NUMBER = 6e23

""" 
Helper functions 
"""
def mass2rad(mass):
   radius = 0.0515*(mass*1000)**(0.393) # Mass in kDa
   return radius #Radius in nm


def rad2mass(radius):
   M = (radius/0.0515)**(1./0.393)/1000.0 #Radius in nm
   return M

def rad2diff(radius):
    viscosity = 0.7e-3 # Pa s
    temperatur = 310.15 # K
    kbT = temperatur*1.38064852e-23
    D = kbT/(6*np.pi*viscosity*radius) #Radius in m
    return D # in m^2/s

class particle:
    def __init__(self, x, S, D, R):
        self.position = x
        self.position0 = x
        self.species = S
        self.diffusion = D
        self.radius = R


class CellList:
    def __init__(self, rc, L, particles):
        self.N_cells = int(L / rc)
        self.cells = [ [[ [] for  col in range(self.N_cells)]
                             for ccol in range(self.N_cells)]
                             for row in range(self.N_cells)]

        self.L = L


    def update(self, particles):
        # Empty cells
        self.cells = [ [[ [] for  col in range(self.N_cells)]
                             for ccol in range(self.N_cells)]
                             for row in range(self.N_cells)]

        cell_ind = [0, 0, 0]
        for i in particles.keys():
            # Place particle in cell lists
            for ind, xij in enumerate(particles[i].position):
                cell_ind[ind] = floor(xij / self.L * self.N_cells)

            self.cells[cell_ind[0]][cell_ind[1]][cell_ind[2]].append(i)

    def get_cells(self, position):
        cell_ind = [0, 0, 0]
        for i, xij in enumerate(position):
            cell_ind[i] = floor(xij / self.L * self.N_cells)
        this_cells = []
        for i in [0, 1]:
            for j in [0, 1]:
                for k in [0, 1]:
                    xind = (cell_ind[0] + i) % self.N_cells
                    yind = (cell_ind[1] + j) % self.N_cells
                    zind = (cell_ind[2] + k) % self.N_cells
                    this_cells.append(self.cells[xind][yind][zind])
                    if not (i == j == k == 0):
                        xind = (cell_ind[0] - i) % self.N_cells
                        yind = (cell_ind[1] - j) % self.N_cells
                        zind = (cell_ind[2] - k) % self.N_cells
                        this_cells.append(self.cells[xind][yind][zind])

        return this_cells

    def check_collisions(self, particles, position, radius):
        intersections = []
        this_cells = self.get_cells(position)
        for this_cell in this_cells:
            for i in this_cell:
                rij = particles[i].position - position
                # Periodic boundaries:
                for k, xij in enumerate(rij):
                    if xij > self.L / 2.:
                        rij[k] = xij - self.L
                    if xij < -self.L / 2.:
                        rij[k] = xij + self.L

                dist = np.sqrt(np.sum((rij) ** 2, axis=0))
                min_dist = particles[i].radius - radius

                if dist < min_dist:
                    intersections.append(i)

        return intersections

    def check_collision(self, particles, position, radius):
        this_cells = self.get_cells(position)
        for this_cell in this_cells:
            for i, p in this_cell.items():
                rij = particles[i].position - position
                # Periodic boundaries:
                for k, xij in enumerate(rij):
                    if xij > self.L / 2:
                        rij[k] = xij - self.L
                    if xij < - self.L / 2:
                        rij[k] = xij + self.L

                dist = np.sqrt(np.sum((rij) ** 2, axis=0))
                min_dist = particles[i].radius - radius

                if dist < min_dist:
                    return True  # collision

        return False  # free


def crowder_free_simulation_method(parameters, phi, seed):
    random.seed(seed)
    np.random.seed(seed)
    # Rescale to mum
    s  = float(1.0e6)
    s2 = float(1.0e12)
    s3 = float(1.0e18)
    # Number of crowders
    R_C = mass2rad(parameters['mu_mass']) * s
    volume_crw = 4.0 / 3.0 * np.pi * R_C ** 3
    N_crw = round(parameters['volume'] / 1000.0 * s3 * phi / volume_crw)
    # Rescaled volume L > mum*3
    V = parameters['volume'] / 1000.0 * s3
    # Box  length
    L = float(parameters['volume'] / 1000.0) ** (1.0 / 3.0) * s
    rc = 35e-3
    #
    gamma = 4.0 * np.pi * (parameters['r_A'] + parameters['r_B']) * (parameters['D_A'] + parameters['D_B']) * s3
    rescaled_keff = parameters['k_fwd'] / AVOGADRO_NUMBER / 1000.0 * s3
    effective_k_binding = gamma * rescaled_keff / (gamma - rescaled_keff)
    effective_k_unbinding = parameters['k_fwd'] * parameters['K_eq']

    particles = {}


    # Bad scaled check collisions
    def check_collision(particles, position, radius):
        for i, p in particles.items():
            rij = particles[i].position - position
            # Periodic boundaries:
            for k, xij in enumerate(rij):
                if xij > L / 2.:
                    rij[k] = xij - L
                if xij < -L / 2.:
                    rij[k] = xij + L

            dist = np.sqrt(np.sum((rij) ** 2.0, axis=0))
            min_dist = particles[i].radius - radius

            if dist < min_dist:
                return True  # collision
        return False  # free

    # Add the species in the concentrations
    N_A = parameters['A_0'] * parameters['volume'] * AVOGADRO_NUMBER
    N_B = parameters['B_0'] * parameters['volume'] * AVOGADRO_NUMBER
    N_C = parameters['C_0'] * parameters['volume'] * AVOGADRO_NUMBER

    print("Add Particles")

    # Add particles
    N = 0
    position = rnd.uniform(0.0, L, 3)
    particles[0] = particle(position,
                            'A',
                            parameters['D_A'] * s2,
                            parameters['r_A'] * s)


    while N < N_A:
        position = rnd.uniform(0, L, 3)
        if not check_collision(particles, position, parameters['r_A'] * s):
            particles[max(particles.keys()) + 1] = particle(position,
                                                            'A',
                                                            parameters['D_A'] * s2,
                                                            parameters['r_A'] * s)
            N += 1

    N = 0
    while N < N_B:
        position = rnd.uniform(0, L, 3)
        if not check_collision(particles, position, parameters['r_B'] * s):
            particles[max(particles.keys()) + 1] = particle(position,
                                                            'B',
                                                            parameters['D_B'] * s2,
                                                            parameters['r_B'] * s)
            N += 1
    N = 0
    while N < N_A:
        position = rnd.uniform(0, L, 3)
        if not check_collision(particles, position, parameters['r_C'] * s):
            particles[max(particles.keys()) + 1] = particle(position,
                                                            'C',
                                                            parameters['D_C'] * s2,
                                                            parameters['r_C'] * s)
            N += 1
    # Recalcualtions
    A = N_crw * R_C / V
    B = 4.0 * np.pi * N_crw * R_C ** 2.0 / V
    C = N_crw * R_C ** 2.0 / V

    def log_p(r):
        return log(1.0 - phi) \
               - B * r / (1.0 - phi) \
               - np.pi * 4.0 * A * r ** 2.0 / (1.0 - phi) \
               - B ** 2.0 * r ** 2.0 / (2.0 * (1.0 - phi) ** 2.0) \
               - 4.0 * np.pi / 3.0 * (N_crw / (V * (1.0 - phi)) + B ** 2.0 * C / (3.0 * (1.0 - phi) ** 3.0) + A * B / (
                    1.0 - phi) ** 2.0) * r ** 3.0

    # P(legal for C)
    p_C = exp(log_p(parameters['r_C'] * s))
    # P(legal for B)
    p_A = exp(log_p(parameters['r_A'] * s))
    # P(legal for A)
    p_B = exp(log_p(parameters['r_B'] * s))

    def sample_spherical(ndim=3):
        vec = np.random.randn(ndim)
        vec /= np.linalg.norm(vec, axis=0)
        return vec

    print('Init Results')

    dt_log = parameters['t_max'] / 10000.0
    t_log = dt_log

    class Result:
        def __init__(self, t, A, B, C):
            self.time = [t, ]
            self.species = {
                'A': [A, ],
                'B': [B, ],
                'C': [C, ],
            }

    result = Result(0, N_A, N_B, N_C)

    print('Init cell list ')


    print(L)
    print(rc)
    # Simulation Algorithm
    cell_list = CellList(rc, L, particles)
    cell_list.update(particles)

    print("Simulation set up")
    t = 0
    while t < parameters['t_max']:
        N = len(particles)

        i, p = random.choice(list(particles.items()))

        # Propagate particles
        dx = rnd.normal(0, sqrt(2.0 * particles[i].diffusion * parameters['dt']))
        particles[i].position += dx

        # Periodic boundaries:
        for k, xij in enumerate(particles[i].position):
            if xij >= L:
                particles[i].position[k] = xij - L
            if xij < 0:
                particles[i].position[k] = xij + L

        dx_dist = np.sqrt(np.sum((dx) ** 2.0, axis=0))

        if rnd.random() < np.pi * N_crw * (R_C + particles[i].radius) * dx_dist / V:
            particles[i].position = particles[i].position0
        # Check for an overlap
        intersections = cell_list.check_collisions(particles,
                                                   particles[i].position,
                                                   particles[i].radius)
        if len(intersections) == 1:
            j = intersections[0]
            try:
                rij = particles[i].position - particles[j].position

                # Periodic boundaries:
                for k, xij in enumerate(rij):
                    if xij > L / 2.0:
                        rij[k] = xij - L
                    if xij < -L / 2.0:
                        rij[k] = xij + L

                dist = np.sqrt(np.sum((rij) ** 2.0, axis=0))
                min_dist = particles[i].radius - particles[j].radius
                if dist < min_dist:
                    if (particles[i].species == 'A' and particles[j].species == 'B') or \
                            (particles[i].species == 'B' and particles[j].species == 'A'):

                        # Reaction monte carlo for product C
                        if rnd.random() < effective_k_binding * parameters['dt']:
                            # P(legal of C)
                            p = p_C

                            if np.random() < 1.0 - p:
                                # Put back
                                particles[i].position = particles[i].position0
                            else:
                                # Check for intersection
                                position = (particles[i].position + particles[j].position) / 2.0
                                if not cell_list.check_collision(particles,
                                                                 position,
                                                                 parameters['r_C'] * s):
                                    # Reaction takes plase
                                    del particles[i]
                                    del particles[j]
                                    particles[max(particles.keys()) + 1] = particle(position,
                                                                                    'C',
                                                                                    parameters['D_C'] * s2,
                                                                                    parameters['r_C'] * s)
                                else:
                                    # Put back
                                    particles[i].position = particles[i].position0

                        else:
                            # Put back
                            particles[i].position = particles[i].position0
                    else:
                        # Put back
                        particles[i].position = particles[i].position0

            except KeyError:
                continue
        else:
            # Put back
            particles[i].position = particles[i].position0

        # Update cell list
        cell_list.update(particles)

        i = min(particles.keys())
        N = max(particles.keys()) + 1
        first_order = False
        while i < N:
            # First order reaction
            try:
                if particles[i].species == 'C':
                    if rnd.random() < effective_k_unbinding * parameters['dt'] / N:
                        p = p_A * p_B

                        if rnd.random() < 1.0 - p:
                            # Reject the reactions
                            pass
                        else:

                            dist = (parameters['r_A'] + parameters['r_B']) * s
                            position1 = particles[i].position + sample_spherical() * dist / 2.0
                            position2 = particles[i].position - sample_spherical() * dist / 2.0
                            # Check intersection
                            if not check_collision(particles, position1, parameters['r_A'] * s) and \
                                    not check_collision(particles, position2, parameters['r_B'] * s):

                                # Periodic boundaries:
                                for k, xij in enumerate(position1):
                                    if xij >= L:
                                        position1[k] = xij - L
                                    if xij < 0:
                                        position1[k] = xij + L

                                # Periodic boundaries:
                                for k, xij in enumerate(position2):
                                    if xij >= L:
                                        position2[k] = xij - L
                                    if xij < 0:
                                        position2[k] = xij + L

                                # Place particles
                                del particles[i]
                                particles[max(particles.keys()) + 1] = particle(position1,
                                                                                'A',
                                                                                parameters['D_A'] * s2,
                                                                                parameters['r_A'] * s)
                                particles[max(particles.keys()) + 1] = particle(position2,
                                                                                'B',
                                                                                parameters['D_B'] * s2,
                                                                                parameters['r_B'] * s)
                                N = max(particles.keys()) + 1
                                first_order = True
            except KeyError:
                i += 1
                continue

            # Next particle
            i += 1

        if first_order:
            # Update cell list again
            cell_list.update(particles)

        t += parameters['dt'] / N
        # print('{}'.format(t))

        if t_log < t:
            result.time.append(t)
            result.species['A'].append(sum([1 for p in particles.values() if p.species == 'A']))
            result.species['B'].append(sum([1 for p in particles.values() if p.species == 'B']))
            result.species['C'].append(sum([1 for p in particles.values() if p.species == 'C']))
            t_log += dt_log

            print('{}'.format(t))
    return result

