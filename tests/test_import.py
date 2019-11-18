""" Tests that the plugin can be loaded and initialized without failure."""

import PDB_plugin as plugin
import pymol

import pytest
import sys


@pytest.fixture(autouse=True)
def InitializePymol():
    # --- setup ---
    # Finish launching PyMOL without GUI.
    if sys.version_info[0] == 2:
        # Under pytest this apparently works only for python2. But it appears
        # that PyMOL under python3 is happy enough without this call.
        pymol.finish_launching(['pymol', '-cqk'])
    # Temporarily turn on maximal log level to also test all logging statements.
    pref_loglevel = 'PDB_PLUGIN_LOGLEVEL'
    orig_loglevel = pymol.plugins.pref_get(pref_loglevel, None)
    pymol.plugins.pref_set(pref_loglevel, 'DEBUG')

    yield  # each test runs here

    # --- teardown ---
    if orig_loglevel is not None:
        pymol.plugins.pref_set(pref_loglevel, orig_loglevel)


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
        '',
        'invalid_pdb_id',
        '3mxw',  # valid PDB ID
        pytest.param('mmCIF_file=tests/data/invalid_file_name',
                     marks=pytest.mark.xfail(raises=pymol.CmdException,
                                             reason='fails to open file')),
        'mmCIF_file=tests/data/3mxw.cif',  # valid file name
        ['3mxw', 'mmCIF_file=tests/data/3mxw.cif'],  # both
    ])
def test_run_main(argv):
    # Load and initialize the module and attempt main().
    # We just want to make sure nothing is crashing at this point.
    # Allow argv singletons for convenience; now wrap them.
    if not isinstance(argv, list):
        argv = [argv]
    plugin.main(argv)
    plugin.count_chains()

    assert True
