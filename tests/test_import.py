""" Tests that the plugin can be loaded and initialized without failure."""

import PDB_plugin as plugin
import pymol

import importlib
import json
import logging
import pytest
import socket
import sys

# ----- Test Fixtures -----

PREF_LOGLEVEL = 'PDB_PLUGIN_LOGLEVEL'


@pytest.fixture(autouse=True)
def initialize_pymol():
    """Properly initialize and set up PyMOL library."""
    # --- setup ---
    # Finish launching PyMOL without GUI.
    if sys.version_info[0] == 2:
        # Under pytest this apparently works only for python2. But it appears
        # that PyMOL under python3 is happy enough without this call.
        pymol.finish_launching(['pymol', '-cqk'])
    # Temporarily turn on maximal log level to also test all logging statements.
    orig_loglevel = pymol.plugins.pref_get(PREF_LOGLEVEL, None)
    pymol.plugins.pref_set(PREF_LOGLEVEL, 'DEBUG')
    # Also set log level directly in logging library for unit tests.
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    pymol.cmd.reinitialize()

    yield  # each test runs here

    # --- teardown ---
    if orig_loglevel is not None:
        pymol.plugins.pref_set(PREF_LOGLEVEL, orig_loglevel)


url_data_cache = {}


@pytest.fixture(autouse=True)
def cache_url_data_for_plugin(monkeypatch):
    """Cache data returned from URLs fetched over the network.

    Wrap URL fetcher cache around code so we don't hit the real PDB with
    repeated requests from this test.
    """
    # --- setup ---
    orig_get_data = getattr(plugin.pdb._fetcher, 'get_data')

    def cached_get_data(url, description):
        if url in url_data_cache:
            return url_data_cache[url]
        data = orig_get_data(url, description)
        url_data_cache[url] = data
        return data

    monkeypatch.setattr(plugin.pdb._fetcher, 'get_data', cached_get_data)

    yield  # each test runs here

    # --- teardown ---
    pass  # nothing to do


# ----- Unit Tests -----


def test_pdb_fetcher_with_requests(requests_mock):
    """Tests the PDB json data fetcher class using requests library."""

    fetcher = plugin.PdbFetcher()

    # Test a response with status code 200 (ie. OK)
    good_response = {'this': ['is', 1, {'good': 'response', 'oh': 'yeah'}]}
    requests_mock.get('http://testpdb/good', json=good_response)
    data = fetcher.get_data('http://testpdb/good', 'good data')
    assert data == good_response

    # Test a response with status code indicating error (e.g. 404 Not Found)
    for status_code in (403, 404, 500):
        requests_mock.get('http://testpdb/bad',
                          status_code=status_code,
                          json='not found')
        fetcher = plugin.PdbFetcher()
        data = fetcher.get_data('http://testpdb/bad', 'no data')
        # We want to see an empty dictionary, nothing else.
        assert isinstance(data, dict)
        assert data == {}


class ImportModuleMock(object):
    """Exclude specified modules from loading in importlib.import_module."""

    def __init__(self, exclude):
        """For excluding all modules specify exclude='*'"""
        self._exclude = exclude
        self._orig_import_module = getattr(importlib, 'import_module')

    # PDB fetcher uses 'requests' library by default. We'll prevent that here
    # such that it falls back to urllib.
    def import_module(self, name):
        if self._exclude == '*' or name in self._exclude:
            raise ImportError('test prevents usage of %s library' % name)
        return self._orig_import_module(name)


def test_pdb_fetcher_with_urllib(monkeypatch):
    """Tests the PDB json data fetcher class using urllib library.

    In the absence of the requests library the code falls  back to urllib2 on
    python2 and urllib.* on python3.
    """

    # PDB fetcher uses 'requests' library by default. We'll prevent that here
    # such that it falls back to urllib.
    mock = ImportModuleMock('requests')
    monkeypatch.setattr(importlib, 'import_module', mock.import_module)
    # We're ready to get a fetcher using urllib.
    fetcher = plugin.PdbFetcher()

    class HTTPResponseMock(object):
        """Returns the stored data json encoded with a read() call."""

        def __init__(self, data):
            self._data = data

        def read(self):
            return json.JSONEncoder().encode(self._data)

    urllib_request = fetcher._modules['urllib.request']
    urllib_error = fetcher._modules['urllib.error']

    # Test a response with status code 200 (ie. OK)
    good_response = {'this': ['is', 2, {'good': 'response', 'oh': 'yeah'}]}
    # Replace urlopen with mock.
    monkeypatch.setattr(urllib_request, 'urlopen',
                        lambda *arg: HTTPResponseMock(good_response))
    data = fetcher.get_data('http://testpdb/good',
                            'good data',
                            sleep_min=0,
                            sleep_max=0)
    assert data == good_response

    def _raise(ex):
        raise ex

    # Test a response with various errors (e.g. status code 404 Not Found).
    for error in (
            urllib_error.HTTPError(None, 403, None, None, None),
            urllib_error.HTTPError(None, 404, None, None, None),
            urllib_error.HTTPError(None, 500, None, None, None),
            urllib_error.URLError(None),
            socket.timeout(),
    ):
        monkeypatch.setattr(urllib_request, 'urlopen',
                            lambda *arg: _raise(error))
        data = fetcher.get_data('http://testpdb/bad',
                                'no data',
                                sleep_min=0,
                                sleep_max=0)
        # We want to see an empty dictionary, nothing else.
        assert isinstance(data, dict)
        assert data == {}


def test_pdb_fetcher_missing_libraries(monkeypatch):
    """Tests the PDB json data fetcher with missing libraries."""

    # Prevent all libraries from getting loaded.
    mock = ImportModuleMock('*')
    monkeypatch.setattr(importlib, 'import_module', mock.import_module)
    # We're ready to get a fetcher using urllib.
    try:
        plugin.PdbFetcher()
    except Exception as e:
        assert 'missing python libraries' in str(e).split('\n')[0]
    else:
        raise Exception('No missing library Exception from PdbFetcher')


def test_pdb_autocomplete(capsys):
    """Tests the PDB ID autocomplete class."""

    # --- Get helpful explanation for valid partial key.
    autocompleter = plugin.PdbIdAutocomplete(get_summary=lambda key: None)
    for i in range(4):
        completion = autocompleter.interpret('1' * i)
        assert completion is not None
        # Completion should be a list > 1 elements to indicate that there are
        # still various options.
        assert isinstance(completion, list)
        assert len(completion) > 1

    # --- Get None for invalid partial key.
    autocompleter = plugin.PdbIdAutocomplete(
        get_summary=lambda key: {key: [{
            'title': 'Title for ' + key
        }]})
    for i in range(1, 4):  # 0-length is covered above
        completion = autocompleter.interpret('a' * i)
        assert completion is None

    # --- Get None for bad key.
    autocompleter = plugin.PdbIdAutocomplete(get_summary=lambda key: None)
    # check cache miss and hit
    for test in ['cache miss', 'cache hit']:
        completion = autocompleter.interpret('bad_key')
        assert completion is None

    # --- Get != None for good key. Also check that returned key is lower()ed.
    autocompleter = plugin.PdbIdAutocomplete(
        get_summary=lambda key: {key: [{
            'title': 'Title for ' + key
        }]})
    key = 'Good_Key'
    assert key != key.lower()  # check for testcase internal bug
    # check cache miss and hit
    for test in ['cache miss', 'cache hit']:
        completion = autocompleter.interpret(key)
        assert completion == key.lower()
        captured = capsys.readouterr()
        assert 'Title for %s' % key.lower() in captured.out

    # --- Check cache purge path.
    autocompleter = plugin.PdbIdAutocomplete(
        get_summary=lambda key: {key: [{
            'title': 'Title for ' + key
        }]})
    # fill up cache with new keys
    for i in range(plugin.PdbIdAutocomplete.MAX_CACHE_SIZE):
        key = 'key%d' % i
        completion = autocompleter.interpret(key)
        assert completion == key
        captured = capsys.readouterr()
        assert 'Title for %s' % key in captured.out
    # exceed cache size; will execute cache purge code path
    key = 'purge_key'
    completion = autocompleter.interpret(key)
    assert completion == key
    captured = capsys.readouterr()
    assert 'Title for %s' % key in captured.out

    # --- Get None when invalid summary format is returned.
    autocompleter = plugin.PdbIdAutocomplete(
        get_summary=lambda key: {key: 'Title for ' + key})
    key = 'good_key'
    completion = autocompleter.interpret(key)
    assert completion is None

    # --- Check built-in access of summary data via PDB lookup.
    autocompleter = plugin.PdbIdAutocomplete()
    key = '3mzw'
    completion = autocompleter.interpret(key)
    assert completion == key
    captured = capsys.readouterr()
    assert 'HER2' in captured.out


def test_initialize():
    """Tests initialization fuction."""
    logger = logging.getLogger()

    # When unset (first time use of pluging), set preference to WARNING.
    pymol.plugins.pref_set(PREF_LOGLEVEL, None)
    plugin.initialize()
    loglevel = pymol.plugins.pref_get(PREF_LOGLEVEL, None)
    assert loglevel == "WARNING"
    assert logger.getEffectiveLevel() == logging.WARNING

    # When set to something valid, keep preference as is.
    pymol.plugins.pref_set(PREF_LOGLEVEL, 'INFO')
    plugin.initialize()
    loglevel = pymol.plugins.pref_get(PREF_LOGLEVEL, None)
    assert loglevel == 'INFO'
    assert logger.getEffectiveLevel() == logging.INFO

    # When set to something invalid, keep preference as is but set actual
    # loglevel to WARNING.
    pymol.plugins.pref_set(PREF_LOGLEVEL, 'BOgus')
    plugin.initialize()
    loglevel = pymol.plugins.pref_get(PREF_LOGLEVEL, None)
    assert loglevel == 'BOgus'
    assert logger.getEffectiveLevel() == logging.WARNING


# ----- Integration Tests -----


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
    """Load and initialize the module and attempt calling main().
    We just want to make sure nothing is crashing at this point.
    """
    # Allow argv singletons for convenience; now wrap them.
    if not isinstance(argv, list):
        argv = [argv]
    plugin.main(argv)
    assert True


def test_commandline_api():
    """Load and initialize the module and attempt calling commandline API.

    Try calling each registered API function to make sure calling them doesn't
    break. We pass the pdbid to each function that is given as an example in
    it's help message. This has the added benefit of verifying that the help
    message example still runs. We do a cursory check on the result of
    count_chains, but otherwise correct functionality is not seriously checked
    here.
    """

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


def test_gui_api(monkeypatch):
    """Load and initialize the module and attempt calling GUI API.

    Try calling each GUI API function registered to a menu item to make sure
    calling them doesn't break. We do a cursory check on the result of
    count_chains, but otherwise correct functionality is not seriously checked
    here.
    """

    gui = plugin.PdbeGui()

    # The GUI normally puts up a dialog to ask for the pdbid. We override that
    # here and provide the dialog result programmatically instead.
    def patch_pdbid_input(pdbid, ok_pressed=True):
        if sys.version_info[0] == 3:
            monkeypatch.setattr(gui._qt.QtWidgets.QInputDialog, 'getText',
                                lambda *arg: (pdbid, ok_pressed))
        else:
            # Above doesn't work under python2, complaining about:
            #     TypeError: unbound method <lambda>() must be called with
            #                QInputDialog instance as first argument (got
            #                NoneType instance instead)
            # So we're patching at a higher level instead.
            monkeypatch.setattr(gui, '_get_pdbid', lambda *arg: pdbid
                                if ok_pressed else None)

    patch_pdbid_input('3l2p')
    gui.analyze_molecules()
    assert plugin.count_chains() == 4

    pymol.cmd.reinitialize()
    patch_pdbid_input('3b43')
    gui.analyze_domains()
    assert plugin.count_chains() == 1

    pymol.cmd.reinitialize()
    patch_pdbid_input('2gc2')
    gui.analyze_validation()
    assert plugin.count_chains() == 2

    pymol.cmd.reinitialize()
    patch_pdbid_input('5j96')
    gui.analyze_assemblies()
    assert plugin.count_chains() == 3

    pymol.cmd.reinitialize()
    patch_pdbid_input('3mzw')
    gui.analyze_all()
    assert plugin.count_chains() == 2

    # Test code doesn't crash when user cancels pdbid input, inputs nothing, or
    # inputs an invalid pdbid.
    pymol.cmd.reinitialize()
    for pdbid, ok_pressed in [('', False), ('', True), ('bogus', False),
                              ('bogus', True)]:
        print('user input', pdbid, ok_pressed)
        patch_pdbid_input(pdbid, ok_pressed)
        gui.analyze_molecules()
        gui.analyze_domains()
        gui.analyze_validation()
        gui.analyze_assemblies()
        gui.analyze_all()


def test_display_types():
    """Tests PDB entries causing use of different display types."""

    # has many polymers -> display type ribbin
    plugin.PDB_Analysis_Molecules('3jcd')
    assert plugin.count_chains() == 57

    # has ca_p_only molecule -> display type ribbon
    pymol.cmd.reinitialize()
    plugin.PDB_Analysis_Molecules('1a1q')
    assert plugin.count_chains() == 3

    # has small polypeptide -> display type sticks
    pymol.cmd.reinitialize()
    plugin.PDB_Analysis_Molecules('6a5j')
    assert plugin.count_chains() == 1

    # has small nucleotide -> display type sticks
    pymol.cmd.reinitialize()
    plugin.PDB_Analysis_Molecules('1b2m')
    assert plugin.count_chains() == 5


def test_gui_missing_libraries(monkeypatch):
    """Tests the GUI with missing libraries."""

    # Prevent all libraries from getting loaded.
    mock = ImportModuleMock('*')
    monkeypatch.setattr(importlib, 'import_module', mock.import_module)
    # We're ready to get a GUI using urllib.
    try:
        plugin.PdbeGui()
    except Exception as e:
        assert 'missing python libraries' in str(e).split('\n')[0]
    else:
        raise Exception('No missing library Exception from PdbFetcher')


def test_object_atom_count():
    """More in-depth test of 'domains' analysis result.

    For one data set, do a more in-depth comparison of generated objects and
    their atom counts.
    """
    plugin.PDB_Analysis_Domains('3mzw')
    assert plugin.count_chains() == 2

    expected_objects = {
        # name: atom_count
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

    # Check all objets' atom counts.
    for object_name in object_names:
        atom_count = pymol.cmd.count_atoms(object_name)
        assert expected_objects[object_name] == atom_count
