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
    pymol.cmd.reinitialize()

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
    assert True


def test_exported_api():
    # Try calling each registered API function to make sure calling them doesn't
    # break. We pass the pdbid to each function that is given as an example in
    # it's help message. This has the added benefit of verifying that the
    # help message example still runs. We do a cursory check on the result of
    # count_chains, but otherwise correct functionality is not seriously checked
    # here.

    plugin.PDB_Analysis_Molecules('3l2p')
    assert plugin.count_chains() == 4

    pymol.cmd.reinitialize()
    plugin.PDB_Analysis_Domains('3b43')
    assert plugin.count_chains() == 1

    pymol.cmd.reinitialize()
    plugin.PDB_Analysis_Validation('2gc2')
    assert plugin.count_chains() == 2

    pymol.cmd.reinitialize()
    plugin.PDB_Analysis_Assemblies('5j96')
    assert plugin.count_chains() == 3


def test_object_atom_count():
    # For one data set, do a more in-depth comparison of generated objects and
    # their atom counts.
    plugin.PDB_Analysis_Domains('3mzw')
    assert plugin.count_chains() == 2

    expected_objects = {
        # name, atom_count
        '3mzw': 4950,
        'Receptor_tyrosineprotein_kinase_erbB2': 4404,
        'Immunoglobulin_Gbinding_protein_A': 430,
        'NACETYLDGLUCOSAMINE': 56,
        'Pfam_PF02216_': 399,
        'Pfam_PF01030_': 1860,
        'Pfam_PF00757_': 1192,
        'Pfam_PF14843_': 693,
        'CATH_3.80.20.20_3mzwA01': 1244,
        'CATH_3.80.20.20_3mzwA03': 1316,
        'CATH_1.20.5.420_3mzwB00': 430,
        'CATH_2.10.220.10_3mzwA02': 817,
        'CATH_2.10.220.10_3mzwA04': 663,
    }

    # Check object names.
    object_names = pymol.cmd.get_object_list()
    assert sorted(expected_objects.keys()) == sorted(object_names)

    # Check atom counts.
    for object_name in object_names:
        atom_count = pymol.cmd.count_atoms(object_name)
        assert expected_objects[object_name] == atom_count
