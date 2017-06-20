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

The topology and docking stages are executed on the local machine and require the software to be available and the
path to the executables defined in the eTOXlie/data/settings.json configuration file.

GROMACS molecular dynamics calculations can be offloaded to an external cluster using an SSH connection.
Configuration of the SSH connection, cluster queueing system and external GROMACS executables should be defined
in the eTOXlie/data/settings.json configuration file.

The eTOXlie installer will try to auto detect prerequisites but a manual check is recommended.


### Installation

The eTOX ALLIES Python code uses a number of third-party Python packages. The easiest way to install eTOX ALLIES
and respect Python dependencies is to use the Python virtual environment. This will install all Python dependencies 
in a Python virtual environment from the Python package repository (pip).

To install using the virtual environment run the eTOX ALLIES installer under bash as:

    >> ./installer.sh -s

The second option is to install eTOX ALLIES without the virtual environment in wich case you are responsible for
installing third-party Python packages yourself. The required Python dependencies are listed in the Pipfile.

To install without the use of the virtual environment run the eTOX ALLIES installer under bash as:

    >> ./installer.sh -s -n

**NOTE**

Please note that the packages installed via pip need the python development headers installed (python-dev)


### Configuration

Application configuration is defined in the eTOXlie/data/settings.json file. Use this file to define eTOXlie application
settings, paths to external dependancies and connection to external clusters for GROMACS calculations.

### Documentation

Documentation on the use of the eTOX ALLIES web based graphical user interface is included in the package as the
eTOX\_ALLIES\_GUI.docx Word document


### Citing eTOX ALLIES

If you have been using eTOX ALLIES please cite it using the following reference:

<cite>L. Capoferri et al. "eTOX ALLIES: an Automated pipeLine for Linear Interaction Energy-based Simulations",
(2017) J. Cheminf.</cite>