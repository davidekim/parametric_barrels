# parametric_barrels

![header.png](./header.png)

## Description
This repo includes utility scripts to generate and evaluate parametrically guided beta barrel protein backbone structures.

## Reference
D.E. Kim, J.L. Watson, D. Juergens, S. Majumder, R. Sonigra, S.R. Gerben, A. Kang, A.K. Bera, X. Li, & D. Baker, Parametrically guided design of beta barrels and transmembrane nanopores using deep learning, Proc. Natl. Acad. Sci. U.S.A. 122 (38) e2425459122, https://doi.org/10.1073/pnas.2425459122 (2025).

## Installation
You can clone this repo into a preferred destination directory by going to that directory and then running:

`git clone https://github.com/davidekim/parametric_barrels.git`

## Usage
barrels.py is the main script that generates parameter defined barrel cylinders used as input for RF partial diffusion and RFJoint2 inpainting beta barrel structure refinement.

`python ./barrels.py --n 6 --S 10 --nres 10`

check_barrels.py is a utility script that evaluates beta barrel outputs from RF partial diffusion and RFJoint2 inpainting to determine the extent of beta barrel formation.

`python ./check_barrels.py <pdblist> <cylinders dir>`

### Dependencies
PyRosetta [https://www.pyrosetta.org](https://www.pyrosetta.org)

BBQ [https://biocomp.chem.uw.edu.pl/tools/bbq](https://biocomp.chem.uw.edu.pl/tools/bbq)

silent_tools [https://github.com/bcov77/silent_tools](https://github.com/bcov77/silent_tools)

## Support
Contact David Kim (dekim@uw.edu) for any questions.


