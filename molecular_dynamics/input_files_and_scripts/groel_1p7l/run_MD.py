#!/usr/bin/env python3
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
from sys import stdout, exit, stderr
import os, time, traceback
import parmed as pmd
import numpy as np

############################################
    
# energy decomposition 
def forcegroupify(system):
    forcegroups = {}
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        force.setForceGroup(i)
        f = str(type(force))
        s = f.split('\'')
        f = s[1]
        s = f.split('.')
        f = s[-1]
        forcegroups[i] = f
    return forcegroups

def getEnergyDecomposition(handle, context, forcegroups):
    energies = {}
    for i, f in forcegroups.items():
        try:
            states = context.getState(getEnergy=True, groups={i})
        except ValueError as e:
            print(str(e))
            energies[i] = Quantity(np.nan, kilocalories/mole)
        else:
            energies[i] = states.getPotentialEnergy()
    results = energies
    handle.write('    Potential Energy:\n')
    for idd in energies.keys():
        handle.write('      %s: %.4f kcal/mol\n'%(forcegroups[idd], energies[idd].value_in_unit(kilocalories/mole))) 
    return results

# remove bond constraints of LIG atoms
def rm_cons_LIG(system, psf_pmd, forcefield, templete_map):
    system_new = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
                 nonbondedCutoff=2.0*nanometer, switchDistance=1.8*nanometer, 
                 constraints=None, removeCMMotion=False, ignoreExternalBonds=True, 
                 residueTemplates=templete_map)
    bond_force = system_new.getForce(0)
    bond_parameter_list = [bond_force.getBondParameters(i) for i in range(bond_force.getNumBonds())]
    tag = 0
    while tag == 0 and system.getNumConstraints() != 0:
        for i in range(system.getNumConstraints()):
            con_i = system.getConstraintParameters(i)[0]
            con_j = system.getConstraintParameters(i)[1]
            segid_i = psf_pmd.atoms[con_i].residue.segid
            segid_j = psf_pmd.atoms[con_j].residue.segid
            if segid_i == 'LIG' and segid_j == 'LIG':
                system.removeConstraint(i)
                # print('Constraint %d is removed, range is %d'%(i, system.getNumConstraints()))
                for bp in bond_parameter_list:
                    if (con_i == bp[0] and con_j == bp[1]) or (con_i == bp[1] and con_j == bp[0]):
                        system.getForce(0).addBond(*bp)
                        break
                tag = 0
                break
            else:
                tag = 1
# END remove bond constraints of LIG atoms

############################################

psffile = sys.argv[1]
corfile = sys.argv[2]
prmfile = sys.argv[3]
temp = float(sys.argv[4])
nproc = '4'
outname = sys.argv[5]
rand = int(sys.argv[6])
sim_steps = int(sys.argv[7])

timestep = 0.015*picoseconds
fbsolu = 0.05/picosecond
temp = temp*kelvin

psf_pmd = pmd.load_file(psffile)
psf_pmd.coordinates = pmd.load_file(corfile).coordinates
psf = CharmmPsfFile(psffile)
cor = CharmmCrdFile(corfile)
forcefield = ForceField(prmfile)
top = psf.topology
# re-name residues that are changed by openmm
for resid, res in enumerate(top.residues()):
    if res.name != psf_pmd.residues[resid].name:
        res.name = psf_pmd.residues[resid].name
templete_map = {}
for chain in top.chains():
    for res in chain.residues():
        templete_map[res] = res.name
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=2.0*nanometer, 
        constraints=AllBonds, removeCMMotion=False, ignoreExternalBonds=True, 
        residueTemplates=templete_map)
custom_nb_force = system.getForce(4)
custom_nb_force.setUseSwitchingFunction(True)
custom_nb_force.setSwitchingDistance(1.8*nanometer)

# Spherical restraint
k = 1.0*kilocalories/mole/angstroms**2
R0 = 160 * angstrom
center_xyz = [0, 0, -74] * angstroms
force = CustomExternalForce("k*dR^2; dR=max(R-R0, 0); R=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2);")
force.addGlobalParameter('k', k)
force.addGlobalParameter('R0', R0)
force.addGlobalParameter('x0', center_xyz[0])
force.addGlobalParameter('y0', center_xyz[1])
force.addGlobalParameter('z0', center_xyz[2])
for atom in psf_pmd.atoms:
    if atom.residue.segid == 'H':
        force.addParticle(atom.idx, [])
system.addForce(force)

# Position restraint
k = 0.1*kilocalories/mole/angstroms**2
force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
force.addPerParticleParameter("k")
force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")
for atom in psf_pmd.atoms:
    if atom.residue.segid != 'H':
        force.addParticle(atom.idx, (k, atom.xx*angstrom, atom.xy*angstrom, atom.xz*angstrom))
system.addForce(force)

# RMSD restraint
K_RMSD = 0.1 * kilocalories/mole/angstroms**2
RMSD0 = 0 * angstroms
restrained_atom_indices = []
for atom in psf_pmd.atoms:
    if atom.residue.segid == 'H':
        restrained_atom_indices.append(atom.idx)
rmsd_cv = RMSDForce(cor.positions, restrained_atom_indices)
energy_expression = 'step(dRMSD) * K_RMSD*dRMSD^2; dRMSD = (RMSD-RMSD0);'
energy_expression += 'K_RMSD = %f;' % K_RMSD.value_in_unit_system(md_unit_system)
energy_expression += 'RMSD0 = %f;' % RMSD0.value_in_unit_system(md_unit_system)
restraint_force = CustomCVForce(energy_expression)
restraint_force.addCollectiveVariable('RMSD', rmsd_cv)
system.addForce(restraint_force)

integrator = LangevinIntegrator(temp, fbsolu, timestep)
integrator.setConstraintTolerance(0.00001)
integrator.setRandomNumberSeed(rand)

# prepare simulation
#platform = Platform.getPlatformByName('CPU')
#properties = {'Threads': nproc}

platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
properties["DeviceIndex"] = "0"

forcegroups = forcegroupify(system)

# Attempt of creating the simulation object (sometimes fail due to CUDA environment)
i_attempt = 0
while True:
    try:
        simulation = Simulation(top, system, integrator, platform, properties)
    except Exception as e:
        print('Error occurred at attempt %d...'%(i_attempt+1))
        traceback.print_exc()
        i_attempt += 1
        continue
    else:
        break

simulation.context.setPositions(cor.positions)
simulation.context.setVelocitiesToTemperature(temp)

# append reporters
simulation.reporters = []
simulation.reporters.append(DCDReporter(outname+'.dcd', 10000, append=False))
simulation.reporters.append(pmd.openmm.reporters.RestartReporter(outname+'.ncrst', 10000, netcdf=True))
simulation.reporters.append(StateDataReporter(outname+'.out', 10000, step=True,
    potentialEnergy=True, temperature=True, progress=False, remainingTime=False,
    speed=True, separator='	'))

getEnergyDecomposition(stdout, simulation.context, forcegroups)
#simulation.minimizeEnergy(tolerance=1*kilojoule/mole)
#getEnergyDecomposition(stdout, simulation.context, forcegroups)

# run production simulation
simulation.step(sim_steps)
print('Done!')
