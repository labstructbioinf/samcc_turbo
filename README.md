# SamCC-Turbo
## version 1.0

Software for automatic detection and measurement of coiled coils in PDB structures.

### Requirements:
1. [DSSP](https://github.com/cmbi/dssp/releases) 
2. [Socket](http://coiledcoils.chm.bris.ac.uk/socket/)
3. Python 3.7.x (recommended 3.7.7)
4. Python dependencies as stated in requirements file
5. (optional) Pymol

### Installation:
1. Download and install DSSP and Socket (links above). Important: SamCC-Turbo does not work with DSSP versions 3.x (format incompatibility issues between DSSP and Socket). We recommend DSSP version 2.3.0
2. Install virtualenv
```bash
pip install virtualenv
```
3. Create new virtual environment
```bash
python -m venv samcc_turbo
```
4. Enter directory of your environment and activate it
```bash
cd samcc_turbo
source bin/activate
```
5. Clone Samcc-Turbo from git
```bash
git clone https://github.com/labstructbioinf/samcc_turbo.git
```
6. Install Python dependencies
```bash
pip install -r requirements.txt
```
7. If you want to use provided example notebooks add virtualenv kernel:
```bash
python -m ipykernel install --name=samcc_turbo
``` 
