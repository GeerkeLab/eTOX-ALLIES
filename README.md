# eTOX ALLIES

eTOX ALLIES is an open source framework that allows the automated prediction of ligand-binding free energies 
requiring the ligand structure as only input. eTOX ALLIES is based on the Linear Interaction Energy (LIE) approach,
an efficient end-point free energy method derived from Free Energy Perturbation (FEP) theory. 
Upon submission of a ligand or dataset of compounds, the tool performs the multiple steps required for binding 
free-energy prediction (docking, ligand topology creation, molecular dynamics simulations, data analysis), 
making use of external open source software where necessary. Moreover, functionalities are also available to 
enable and assist the creation and calibration of new models. 
In addition, a web graphical user interface (GUI) has been developed to allow use of free-energy based models to
users that are not an expert in molecular modeling.


### Prerequisites

eTOX ALLES depends on external software for the key stages in the ligand-binding free energy predictions:

* Topology generation: [AMBERTOOLS](http://ambermd.org) and [ACPYPE](https://github.com/t-/acpype)
* Docking: [PLANTS](http://www.tcd.uni-konstanz.de/research/plants.php) or [ParaDocks](https://github.com/cbaldauf/paradocks)
* Molecular Dynamics: [GROMACS](http://www.gromacs.org)

The topology and docking stages are executed on the local machine and require the software to be available in the
system PATH or user defined in the eTOXlie/data/settings.json configuration file.
GROMACS molecular dynamics calculations can be offloaded to an external cluster using an SSH connection.

The eTOXlie installer will try to auto detect prerequisites but a manual check is recommended.


### Installation

Run the eTOX ALLIES installer under bash:

    >> ./installer.sh -s

This will install all Python dependencies in a Python virtual environment from the Python package repository (pip).


### Configuration

Application configuration is defined in the eTOXlie/data/settings.json file. Use this file to define eTOXlie application
settings, paths to external dependancies and connection to external clusters for GROMACS calculations.


### Citing eTOX ALLIES

If you have been using eTOX ALLIES please cite it using the following reference:

<cite>L. Capoferri et al. "eTOX ALLIES: an Automated pipeLine for Linear Interaction Energy-based Simulations",
(2017) J. Cheminf.</cite>