""" Tests that the plugin can be loaded and initialized without failure."""

import PDB_plugin as plugin


def test_initialize_plugin():
    # Initialize the plugin, but don't perform any actions. The initialization
    # can't actually 'addmenuitemqt' since pymol.qt hasn't been initialized, but
    # that's OK, we just want to check that the python code in the plugin
    # initialization doesn't crash.
    plugin.__init_plugin__()
    assert True


def test_run_main_without_arguments():
    # Load and initialize the module and attempt main(). This will receive an
    # error message about missing command line arguments. That's fine, we just
    # want to make sure nothing is crashing at this point.
    argv = []
    plugin.main(argv)
    assert True
