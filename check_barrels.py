#!/usr/bin/env python3
import os, sys
import numpy as np
import pickle
import re
from pathlib import Path
home = str(Path.home())

#
# Reference for parametrically guided beta barrel backbone design:
#   Kim DE. et al. 2024. Parametrically guided design of soluble beta barrels and transmembrane nanopores using deep learning.
#

# Dependencies
# https://github.com/bcov77/silent_tools - silent_tools
# https://www.pyrosetta.org - PyRosetta

import pyrosetta
from pyrosetta import rosetta

# point this to your silent_tools location
sys.path.append( home+'/src/silent_tools' )
import silent_tools

verbose = False

## This script checks for barrel-like topology of parametrically guided barrel designs from the
## number of backbone H-bonds of each adjacent strand pair. For increasing 0.10 incremental cutoffs,
## if the number of H-bonds is at least the fractional cutoff of the strand length for each strand
## pair including the N to C strand pair (closed barrel), the structure passes as a barrel topology.
## As the cutoff increases to 1 (number of H-bonds at least the strand length), the structures contain
## more strands, typically. Structures that only pass at 0.10 contain at least one strand pair with 
## strand length * 0.10 H-bonds.

if len(sys.argv) != 3 and len(sys.argv) != 1:
    print(f'USAGE: {sys.argv[0]} <pdblist - list of barrel design PDBs> <cylinders dir - directory where bb cylinder files exist>')
    print()
    print(f'Outputs <pdblist.silent and pdblist.barrel_status> of barrels')
    exit(1)

def score_pose(p):
    scorefxn_tal_name="beta_nov15"
    scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function(scorefxn_tal_name)
    energymethodoptions=scorefxn.energy_method_options()
    energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(True)
    scorefxn.set_energy_method_options( energymethodoptions );
    return scorefxn(p)

def find_hbonds(p,ssstr):
    global hbs,hbsacc
    # initialize hbs and hbsacc
    for i,aa in enumerate(ssstr):
        iplus = i+1
        if iplus not in hbs:
            hbs[iplus] = []
        if iplus not in hbsacc:
            hbsacc[iplus] = []
    x=score_pose(p)
    p.update_residue_neighbors();
    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
    hbond_set.setup_for_residue_pair_energies(p , False, False);
    for i_hbond in range(hbond_set.nhbonds() ):
        test_hbond=hbond_set.hbond(i_hbond+1)
        if (test_hbond.acc_atm_is_backbone() and test_hbond.don_hatm_is_backbone()):
            accRes=test_hbond.acc_res()
            donRes=test_hbond.don_res()
            atomAname=p.residue(test_hbond.acc_res()).atom_name(test_hbond.acc_atm())
            atomBname=p.residue(test_hbond.don_res()).atom_name(test_hbond.don_hatm())
            if ssstr[int(donRes)-1] == 'E' and  ssstr[int(accRes)-1] == 'E':
                hbsacc[int(accRes)].append(int(donRes))
                hbs[int(donRes)].append(int(accRes))

pdblist = ''
barrelpdbs = []
cylinder_path = './'
if len(sys.argv) == 1:
    import glob
    pdblist = sys.argv[0]
    pdblist = pdblist.split('.py')[0]
    for pdb in glob.glob('n*S*nres*/*.pdb'):
        barrelpdbs.append(pdb)
else:
    pdblist = sys.argv[1]
    with open(pdblist) as pdbl:
        for pdb in pdbl:
            pdb = pdb.strip()
            barrelpdbs.append(pdb)
    cylinder_path = sys.argv[2]

pyrosetta.init('-beta_nov16 -in:file:silent_struct_type binary -out:file:silent_struct_type binary -mute all -out::file::renumber_pdb -in::ignore_unrecognized_res -in::ignore_waters')

outsilent = pdblist+'.silent'

sfd_out = rosetta.core.io.silent.SilentFileData(outsilent, False, False, "binary", rosetta.core.io.silent.SilentFileOptions())

def add2silent( pose, tag, sfd_out ):
    struct = sfd_out.create_SilentStructOP()
    struct.fill_struct( pose, tag )
    sfd_out.add_structure( struct )
    sfd_out.write_silent_struct( struct, outsilent )


def read_span(span):
    [start,stop] = span.split('-')
    start = int(start[1:])
    stop = int(stop)
    return [start,stop]

if not cylinder_path.endswith('/'):
    cylinder_path = cylinder_path + '/'

nf = open(pdblist+'.barrel_status', 'w')
nf.write(f'pdb group n S len min_frac_adj_strand_hb pass_as_barrel cylinder_rmsd sampled_mask'+"\n")
for pdb in barrelpdbs:
    nflinebest = ""
    pdb = pdb.strip()
    cylinder = pdb.split('-bb')[0]+'-bb.pdb' 
    if not os.path.exists(cylinder):
        cylinder = cylinder_path + cylinder.split('/')[-1]
    print(pdb)
    print(cylinder)

    # get sampled mask if inpainted
    trb = pdb.split('.pdb')[0]+'.trb'
    sampled_mask = []
    sampled_mask_str = ""
    if os.path.exists(trb):
        f = open(trb, 'rb')
        trbdat = pickle.load(f)
        sampled_mask = trbdat['sampled_mask'][0].split(',')
        sampled_mask_str = trbdat['sampled_mask'][0]

    hbs = {} 
    hbsacc = {}

    m = re.search('_n([\d]+)_S([\d]+)_nres([\d\.]+)_(.+)_looplen([\d\.]+)_terminilen([\d\.]+)_', pdb)
    if m:
        strands = int(float(m.group(1)))
        shear = int(float(m.group(2)))
        strandlen = int(float(m.group(3)))
        looplen = int(float(m.group(5)))
        terminilen = int(float(m.group(6)))

        if len(sampled_mask) <= 1:
            sampled_mask = [ str(terminilen)+'-'+str(terminilen) ]
            for i in range(1,strands+1):
                span = 'A'+str((strandlen*i-strandlen)+1)+'-'+str(i*strandlen)
                sampled_mask.append(span)
                if i != strands:
                    sampled_mask.append( str(looplen)+'-'+str(looplen) )
            sampled_mask.append(str(terminilen)+'-'+str(terminilen))
            sampled_mask_str = ",".join(sampled_mask)

        print(sampled_mask_str)

        p = pyrosetta.pose_from_file(pdb)

        # get DSSP
        DSSP = pyrosetta.rosetta.core.scoring.dssp.Dssp(p)
        ssstr = DSSP.get_dssp_secstruct()
            
        find_hbonds(p, ssstr)
        best = ""
        for min_frac_adj_strand_hb in [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]:
            pass_as_barrel = True
            hbondcnt = 0

            # try to make sure adjacent strands have enough? hbonds to complete a full barrel structure
            resnumshift = 0
            strand_resnums = []
            for i in range(1,len(sampled_mask)-1,2):
                resnumshift = resnumshift + int(sampled_mask[i-1].split('-')[0])
                [istart,istop] = read_span(sampled_mask[i])
                adjstrandindex = i+2
                jresnumshift = resnumshift + int(sampled_mask[i+1].split('-')[0])
                if i == len(sampled_mask)-2:
                    adjstrandindex = 1
                    jresnumshift = int(sampled_mask[0].split('-')[0])

                [jstart,jstop] = read_span(sampled_mask[adjstrandindex])
                hbondcnt = 0
                for x in range(istart+resnumshift,istop+resnumshift+1):
                    strand_resnums.append(x)
                    for y in range(jstart+jresnumshift,jstop+jresnumshift+1):
                        if y in hbs[x] and x in hbsacc[y]:
                            hbondcnt = hbondcnt + 1
                        if x in hbs[y] and y in hbsacc[x]:
                            hbondcnt = hbondcnt + 1
                if hbondcnt < strandlen*min_frac_adj_strand_hb:
                    pass_as_barrel = False
                rmsd = 0.0    
                if i == len(sampled_mask)-2:
                    # figure out RMSD to cylinder here
                    rmap = pyrosetta.rosetta.std.map_unsigned_long_unsigned_long()
                    for i,r in enumerate(strand_resnums):
                        rmap[i+1] = r
                    c = pyrosetta.pose_from_file(cylinder)
                    rmsd = pyrosetta.rosetta.core.scoring.CA_rmsd( c, p, rmap)
                if verbose:
                    print(f'slen: {strandlen} hbcnt: {hbondcnt} ist: {istart+resnumshift} {istop+resnumshift} jst: {jstart+jresnumshift} {jstop+jresnumshift} ishift: {resnumshift} jshift: {jresnumshift} adjsi: {adjstrandindex} {sampled_mask[i]}')
            barrel_group = f'n{strands}_S{shear}_len{strandlen}' 
            nfline = f'{pdb} {barrel_group} {strands} {shear} {strandlen} {min_frac_adj_strand_hb} {pass_as_barrel} {rmsd} {sampled_mask_str}'+"\n"
            nf.flush()

            if pass_as_barrel:
                best = pdb.split('/')[-1].split('.pdb')[0]
                best = best+f'_{min_frac_adj_strand_hb}_{rmsd:2.2f}'
                nflinebest = nfline
            else:
                nf.write(nfline)
        if len(best) > 0:
            print(best)
            nf.write(nflinebest)
            add2silent( p, best, sfd_out )

nf.close()









