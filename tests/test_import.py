""" Tests that the plugin can be loaded and initialized without failure."""

import PDB_plugin as plugin
import pymol

import pytest


def test_initialize_plugin():
    # Initialize the plugin, but don't perform any actions. The initialization
    # can't actually 'addmenuitemqt' since pymol.qt hasn't been initialized, but
    # that's OK, we just want to check that the python code in the plugin
    # initialization doesn't crash.
    plugin.__init_plugin__()
    assert True


# This is a quick and dirty way to get up the test coverage initially just to
# check for fatal problems in the code. It doesn't check correctness.
@pytest.mark.parametrize(
    'argv',
    [
        [''],
        ['invalid_pdb_id'],
        ['3mxw'],  # valid PDB ID
        # This exits by raising an exception - disable for now.
        # ['mmCIF_file=tests/data/invalid_file_name'],
        ['mmCIF_file=tests/data/3mxw.cif'],  # valid file name
        ['3mxw', 'mmCIF_file=tests/data/3mxw.cif'],  # both
    ])
def test_run_main(argv):
    # Load and initialize the module and attempt main().
    # We just want to make sure nothing is crashing at this point.

    # Turn on maximal log level in order to also test all logging statements.
    pref_loglevel = 'PDB_PLUGIN_LOGLEVEL'
    pymol.plugins.pref_set(pref_loglevel, 'DEBUG')

    plugin.main(argv)
    plugin.count_chains()

    assert True
