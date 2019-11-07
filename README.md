# PDB plugin for PyMOL

This repository is based on [PDBe's](https://pdbe.org/) PDB plugin for
[PyMOL](https://pymol.org/). The
[original documentation](https://pymolwiki.org/index.php/PDB_plugin)
for the plugin is found on the
[PyMOL Wiki](https://pymolwiki.org/index.php/PDB_plugin).

## Versions

The plugin has been updated to newer versions of PyMOL and Python. Here is
an overview of the download URLs for various versions:
* PyMOL 2.x and Python 3:
  * https://github.com/ArturKlauser/pymol-pdb-plugin/raw/master/PDB_plugin.py
* PyMOL 1.x and Python 3:
  * https://github.com/ArturKlauser/pymol-pdb-plugin/raw/v1.x-python3/PDB_plugin.py
  * *The original plugin doesn't work for Python 3 based installations due to
changed locations of tk imports.*
* PyMOL 1.x and Python 2 (original plugin from PDBe):
  * This repository: https://github.com/ArturKlauser/pymol-pdb-plugin/raw/v1.x/PDB_plugin.py
  * Original URL: http://www.ebi.ac.uk/pdbe/pdb-component-library/pymol_plugin/PDB_plugin.py

## Installation
To install the plugin for PyMOL 1.7 or later use the PyMOL menu:
  * *Plugin* -> *Plugin Manager*
  * select the *Install New Plugin* tab
  * fill in the *URL* field with one of the URLs above
  * click *Fetch*
  
For earlier versions of PyMOL install the plugin manually from the above URL.
