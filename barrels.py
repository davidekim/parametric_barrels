#! /usr/bin/env python
#
# Reference for parametrically guided beta barrel backbone design:
#   Kim DE. et al. 2024. Parametrically guided design of soluble beta barrels and transmembrane nanopores using deep learning.
#
# References for cylinder generation code:
#   Naveed H. et al. JACS, 2012. Predicting three-dimensional structures of transmembrane domains of beta-barrel membrane proteins.
#   Dou, J. et al. Nature, 2018. De novo design of a fluorescence-activating β-barrel. 
#
#
# Dependencies
# https://www.pyrosetta.org - PyRosetta
# https://biocomp.chem.uw.edu.pl/tools/bbq - BBQ

from pathlib import Path
home = str(Path.home())

# EDIT THE FOLLOWING TO POINT TO YOUR SPECIFIC PATHS
RUN_BBQ_SH = home+"/src/bioshell.bioinformatics-2.2/run_bbq.sh"
INPAINT_PY = "/software/containers/SE3nv.sif "+home+"/src/proteininpainting/inpaint.py"
RUN_INFERENCE_PY = "/software/containers/RF_diffusion.sif /projects/ml/rf_diffusion/run_inference.py"


##
##
##

import os,sys
from math import sqrt,asin,pi,sin,tan,cos
from scipy.optimize import fsolve
import random
import argparse
from pyrosetta import *
from pyrosetta.rosetta import *
from rosetta.protocols.loops.loop_closure.ccd import *

# Relevant parameter studies
# Cyril F. Reboul, Khalid Mahmood, James C. Whisstock, Michelle A. Dunstone, Predicting
# giant transmembrane β-barrel architecture, Bioinformatics, Volume 28, Issue 10, 15 May 2012,
#  Pages 1299–1302, https://doi.org/10.1093/bioinformatics/bts152

# symmetric membrane
# 3.48 CA-CA intra-strand distance
# 4.83 CA-CA inter-strand

# water-soluble (Chou and Scheraga, 1982)
# 3.53
# 4.87  

# ECOD40 outer membrane meander beta-barrels
# 3.8
# 4.8

parser = argparse.ArgumentParser(description='Program')
parser.add_argument('--n', action="store",type=int, default=8, help="Number of beta strands of beta-barrel")
parser.add_argument('--S', action="store",type=int, default=10, help="Shear number of beta-barrel")
parser.add_argument('--nres', action="store", default=10, type=int,help="number of residues per strand of beta-barrel")
parser.add_argument('--b', action="store",default=4.8, required=False,  type=float,help="CA-CA inter strand distance") 
parser.add_argument('--a', action="store",default=3.8, required=False,  type=float,help="CA-CA intra strand distance") 
parser.add_argument('--terminilen', action="store",default=2, required=False, type=int,help="termini length") 
parser.add_argument('--looplen', action="store",default=4, required=False, type=int,help="loop length") 
parser.add_argument('--inpaint_template_conf', action="store",default=0.9, required=False, type=float,help="Template confidence for inpainting") 
parser.add_argument('--diffusion_partial_T', action="store",default=10, required=False, type=int,help="Partial diffusion partial_T")
parser.add_argument('--dtw', action="store", default=1.0, required=False, type=float,help="tilt angle ratio")
parser.add_argument('--antiparallel', type=bool, default=True, action="store", required=False, help="Anti-parallel")
args = parser.parse_args()

a = args.a #Distance between consective C-alpha atoms within a strand
b = args.b #Distance between C-alpha atoms of adjacent strands which forms hbonds
n = args.n #Number of beta strands of beta-barrel
S = args.S #Shear number of beta-barrel
nres = args.nres #Number of residues per beta strand in beta-barrel
termlen = args.terminilen #Termini length to model (0-termlen for inpainting, termlen for partial diffusion)
looplen = args.looplen #Loop length to model between beta strands (0-looplen for inpainting, looplen for partial diffusion)
template_conf = args.inpaint_template_conf #Template confidence for inpainting input
partial_T = args.diffusion_partial_T #partial_T value for partial diffusion

print(args)

init('-detect_disulf false -out::file::renumber_pdb -in::ignore_unrecognized_res -in::ignore_waters')

outpdb = f'barrelCylinder_n{n}_S{S}_nres{nres}_b{b}_dtw{args.dtw}_looplen{args.looplen}_terminilen{args.terminilen}_cylinder.pdb' 

PI = pi
r=sqrt((n*b)**2+(S*a)**2)/(2*PI)
theta=asin(S*a/(2*PI*r))*args.dtw
def disNextRes(x):
    return sqrt(r**2*(2-2*cos(x))+(r*x/tan(theta))**2)-a
def disNextStrand(x):
    return sqrt(r**2*(2-2*cos(x+2*PI/n))+(r*x/tan(theta))**2)-b
delta_t1 = fsolve(disNextRes,0)[0]
delta_t2 = fsolve(disNextStrand,0)[0]
def dis(x,y):
    s=0
    for i in range(len(x)):s += (x[i]-y[i])**2
    return sqrt(s)
coor = []
for ns in range(1,n+1):
    phi = (ns-1)*2*pi/n
    dt2 = (delta_t2)*(ns-1)
    coor_strand = []
    n0 = 0 #number of residues at the strand
    for j in range(800):
        dt1 = (delta_t1)*j
        dt = dt1+dt2
        x=r*cos(dt+phi)
        y=r*sin(dt+phi)
        z=r*dt/tan(theta)
        if z<0:continue
        n0 += 1
        if n0>nres: break # +4: break #jump out of the loop
        sign = 1 if j%2 == 0 else -1
        coor_strand.append([x,y,z,sign])
    coor.append(coor_strand)
foA = open(outpdb,'w')
resi = 1
atomi = 2
chaini = 0
for i in range(len(coor)):
    coor_strand = coor[i]
    if args.antiparallel and i%2 == 0: coor_strand.reverse() 
    for j in range(len(coor_strand)):
        (x,y,z,sign) = coor_strand[j]
        foA.write("ATOM  %5d  CA  %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %1d\n" % (atomi,"ALA",'A',resi,x,y,z,sign))
        resi += 1	
        atomi += 4
    chaini += 1
foA.close()

# Create backbone atoms and loops for inpainting and/or diffusion input
# place backbone atoms w/ bbq
os.system(f"source {RUN_BBQ_SH} "+outpdb)
bbA = outpdb.split('.pdb')[0]+'-bb.pdb'
prefix = bbA[0:-4]

pose = pose_from_pdb(bbA)
singlechain = Pose()
for i in range(1, len(pose.sequence())+1):
    core.conformation.remove_upper_terminus_type_from_conformation_residue(pose.conformation(), i)    
    core.conformation.remove_lower_terminus_type_from_conformation_residue(pose.conformation(), i)    
    singlechain.append_residue_by_bond(pose.residue(i))

EXT_PHI = -150.0
EXT_PSI = +150.0
EXT_OMG = 180.0
def perturb_bb(pose, pos):
    phi = EXT_PHI + random.randint(-60, 60)
    psi = EXT_PSI + random.randint(-60, 60) 
    omg = EXT_OMG
    pose.set_phi(pos, phi)
    pose.set_psi(pos, psi)
    pose.set_omega(pos, omg)

pose = Pose(singlechain)

# add loops
inserted = 0
loops = protocols.loops.Loops()
for i in range(1,n):  # strands
    previous = i*nres+inserted
    for j in range(looplen):
        insertpos = i*nres+inserted
        res_type = core.chemical.ChemicalManager.get_instance().residue_type_set( 'fa_standard' ).get_representative_type_name1('A')
        residue = core.conformation.ResidueFactory.create_residue(res_type) 
        pose.append_polymer_residue_after_seqpos(residue, insertpos, True)
        inserted += 1
    cutpoint = previous+looplen
    loop = protocols.loops.Loop(previous, cutpoint+1, cutpoint, 0.0, True)
    loops.add_loop(loop)
    protocols.loops.set_single_loop_fold_tree(pose, loop)
    # extend loops before closing
    for l in range(loop.start()+1,loop.stop()):
        core.conformation.idealize_position(l, pose.conformation())
        pose.set_phi(l, EXT_PHI)
        pose.set_psi(l, EXT_PSI)
        pose.set_omega(l, EXT_OMG)
# close loops
lm = protocols.loop_modeler.LoopModeler()
lm.set_loops(loops)
lm.enable_build_stage()
lm.disable_centroid_stage()
lm.disable_fullatom_stage()
lm.apply(pose)

# add extended termini with some phi/psi randomness
for i in range(termlen):
    res_type = core.chemical.ChemicalManager.get_instance().residue_type_set( 'fa_standard' ).get_representative_type_name1('A')
    residue = core.conformation.ResidueFactory.create_residue(res_type)
    core.conformation.idealize_position(1, pose.conformation())
    pose.prepend_polymer_residue_before_seqpos(residue, 1, True)
    core.conformation.idealize_position(1, pose.conformation())
    perturb_bb(pose, 1)
    core.conformation.idealize_position(len(pose.sequence()), pose.conformation())
    perturb_bb(pose, len(pose.sequence()))
    pose.append_polymer_residue_after_seqpos(residue, len(pose.sequence()), True)
    perturb_bb(pose, len(pose.sequence()))

pose.dump_pdb(f'{prefix}_with_loops.pdb')
pose_with_loops_len = len(pose.sequence())

# create inpainting command
contigs = f'0-{termlen},'
end = 0
nresorig = nres
for i in range(1,n+1):
    start = int(nresorig*i-nresorig+1)
    end = int(nresorig*i)
    contigs = contigs+f'A{start}-{end},0-{looplen},'
contigs = contigs[:-1]
contigslist = contigs.split(',')[0:-1]
contigs = ','.join(contigslist)
contigs = contigs + f',0-{termlen}'
subd = f'n{n}_S{S}_nres{nres}_inpainting'
if not os.path.exists(subd):
    os.system(f'mkdir {subd}')
output = bbA.split(".pdb")[0]+f'_looplen_{looplen}_termlen_{termlen}_tmpl_conf_{template_conf}'
inpaint_cmd = f'{INPAINT_PY} --pdb {bbA} --out {subd}/{output} --contigs {contigs} --tmpl_conf {template_conf} --inpaint_seq A1-{end} --num_designs 1'
with open(output+'_inpaint.sh', "w") as inps:
    inps.write(inpaint_cmd+"\n")

# create partial diffusion command
subd = f'n{n}_S{S}_nres{nres}_partial_diffusion'
if not os.path.exists(subd):
    os.system(f'mkdir {subd}')
contigstr = f'{pose_with_loops_len}-{pose_with_loops_len}'
output = bbA.split(".pdb")[0]+f'_looplen_{looplen}_termlen_{termlen}_partial_T_{partial_T}'
partial_diffusion_cmd = f'{RUN_INFERENCE_PY} inference.output_prefix={subd}/{output} inference.input_pdb={prefix}_with_loops.pdb diffuser.partial_T={partial_T} '
partial_diffusion_cmd += f'contigmap.contigs=[\\\'{contigstr}\\\'] inference.num_designs=1 denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5 inference.num_designs=1'
with open(output+'_partial_diffusion.sh', "w") as inps:
    inps.write(partial_diffusion_cmd+"\n")



