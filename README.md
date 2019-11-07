# PDB plugin for PyMOL

This repository is based on [PDBe's](https://pdbe.org/) PDB plugin for
[PyMOL](https://pymol.org/). The
[original documentation](https://pymolwiki.org/index.php/PDB_plugin)
for the plugin is found on the
[PyMOL Wiki](https://pymolwiki.org/index.php/PDB_plugin).
For fixes for Python 3 see below.

## Installation
To install the plugin for PyMOL 1.7 or later use the PyMOL menu:
  * *Plugin* -> *Plugin Manager*
  * select the *Install New Plugin* tab
  * fill in the *URL* field:
    * Original PDBe URL: http://www.ebi.ac.uk/pdbe/pdb-component-library/pymol_plugin/PDB_plugin.py
    * URL to this repository: https://github.com/ArturKlauser/pymol-pdb-plugin/raw/v1.x/PDB_plugin.py
  * click *Fetch*
  
For PyMOL 1.6 install the plugin manually from the above URL.

## Python 3
The original plugin doesn't work for Python 3 based installations due to
changed locations of tk imports. The fixed plugin can be installed from URL 
https://github.com/ArturKlauser/pymol-pdb-plugin/raw/v1.x-python3/PDB_plugin.py
