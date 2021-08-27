# SamCC-Turbo
## version 0.0.2

Software for automatic detection and measurement of coiled coils in PDB structures.

### Requirements
1. [DSSP 2.3.0](https://github.com/cmbi/dssp/releases)
2. [Socket](http://coiledcoils.chm.bris.ac.uk/socket/)
3. Python 3.7.x (recommended 3.7.7)
4. (optional) Pymol. For installation instructions see <https://pymolwiki.org/index.php/Linux_Install>

### Installation
1. Download and install DSSP and Socket (links above). Important: SamCC-Turbo does not work with DSSP versions 3.x (format incompatibility issues between DSSP and Socket). We recommend DSSP version 2.3.0
2. Install virtualenv
```bash
pip install virtualenv
```
3. Create new virtual environment. Please ensure that you are using Python 3. You may need to install python3-venv package. Remember that when creating virtualenv you can specify version of Python to be used with --python argument (we currently recommend Python 3.7.7)
```bash
python -m venv samcc_turbo
```
4. Enter directory of your environment and activate it
```bash
cd samcc_turbo # virtual environment
source bin/activate
```
5. Install SamCC-Turbo from pip
```bash
pip install samcc-turbo
```

### Commands
1. Basic usage (measure PDB, do not save pymol session)
```python
from samcc import SamCC
sc = SamCC(bin_paths={'dssp':'/path/to/dssp', 'socket':'/path/to/socket'})
sc.run_samcc_turbo(pdbpath='/path/to/pdb', outpath='/path/to/results',  save_pse=False) # to save pymol session change save_pse to True
```

2. Commands: to obtain full list of commands for SamCC-Turbo:
```python
from samcc import SamCC
help(SamCC)
```

3. For advanced example see the <Example.ipynb> notebook!

### Do you prefer a Web server?
Please visit <https://lbs.cent.uw.edu.pl/ccdb> and <https://lbs.cent.uw.edu.pl/samcc_turbo>

### Interested in coiled coils?
Please check our papers:
* [DeepCoil-a fast and accurate prediction of coiled-coil domains in protein sequences.](https://www.ncbi.nlm.nih.gov/pubmed/30601942)
* [Variability of the core geometry in parallel coiled-coil bundles.](https://www.ncbi.nlm.nih.gov/pubmed/30042011)
* [The Structure and Topology of α-Helical Coiled Coils.](https://www.ncbi.nlm.nih.gov/pubmed/28101860)
* [Designability landscape reveals sequence features that define axial helix rotation in four-helical homo-oligomeric antiparallel coiled-coil structures.](https://www.ncbi.nlm.nih.gov/pubmed/25278129)

### Questions?
If you encounter any problems with the installation or running SamCC-Turbo, please do not hesitate to contact us at <s.dunin-horkawicz@cent.uw.edu.pl>
