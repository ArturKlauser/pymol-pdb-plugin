""" Tests that the plugin can be loaded and initialized without failure."""

import PDB_plugin as plugin
import pymol

import importlib
import json
import logging
import os
import pytest
import socket
import sys
try:
    import urllib.parse as url_parse
except ImportError:
    import urllib as url_parse

# ----- Test Fixtures -----

PREF_LOGLEVEL = 'PDB_PLUGIN_LOGLEVEL'
WEBCACHE_PATH = 'tests/data/webcache'


@pytest.fixture(autouse=True)
def initialize_pymol():
    """Properly initialize and set up PyMOL library."""
    # --- setup ---
    # Finish launching PyMOL without GUI.
    if sys.version_info[0] == 2:
        # Under pytest this apparently works only for python2. But it appears
        # that PyMOL under python3 is happy enough without this call.
        pymol.finish_launching(['pymol', '-cqk'])
    pymol.cmd.reinitialize()
    pymol.plugins.initialize(pmgapp=-2)  # No Autoloading
    # Don't save preference changes into rc file.
    pymol.plugins.pref_set('instantsave', False)
    # Temporarily turn on maximal log level to also test all logging statements.
    pymol.plugins.pref_set(PREF_LOGLEVEL, 'DEBUG')
    # Also set log level directly in logging library for unit tests.
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    yield  # each test runs here

    # --- teardown ---
    pass  # nothing to do


class WebCache(object):
    """Cache web data accessed by the PDB API URLs."""

    def __init__(self, fetcher):
        """Uses 'fetcher' to get data on cache miss."""
        self._fetcher = fetcher if fetcher else self._get_fetch_exception
        self._cache = {}
        self._num_accesses = 0
        self._num_hits = 0

    @staticmethod
    def _get_fetch_exception(url, description):
        raise Exception(
            'Missing in webcache: %s (%s)\nFetching from web is disallowed.' %
            (description, url))

    def access(self, url, description):
        """Returns the data for the given URL."""
        self._num_accesses += 1
        if url in self._cache:
            self._num_hits += 1
            return self._cache[url]
        data = self._fetcher(url, description)
        self._cache[url] = data
        return data

    def save_to_dir(self, path):
        """Saves cache data, one file per key, located in dir path."""
        print('Saving %d URLs to webcache %s.' % (len(self._cache), path))
        if not os.access(path, os.O_DIRECTORY):
            os.mkdir(path)
        for url, data in self._cache.items():
            url = url_parse.quote_plus(url, safe='')  # must quote '/'
            with open(os.path.join(path, url), 'w') as file:
                json.dump(data, file, indent=2, sort_keys=True)

    def load_from_dir(self, path):
        """Loads cache data, one file per key, located in dir path."""
        print('loading webcache from', path)
        for url in os.listdir(path):
            with open(os.path.join(path, url), 'r') as file:
                url = url_parse.unquote_plus(url)
                self._cache[url] = json.load(file)

    def report_stats(self):
        print('\nWebCache Stats: %d accesses, %d hits, %d misses; '
              '%3.0f%% hit rate' %
              (self._num_accesses, self._num_hits,
               self._num_accesses - self._num_hits, 100.0 * (self._num_hits) /
               self._num_accesses if self._num_accesses else 0))


@pytest.fixture(scope='session')
def web_cache(pytestconfig):
    # --- setup ---
    if pytestconfig.option.webcache_fetch:
        fetcher = getattr(plugin.pdb._fetcher, 'get_data')
    else:
        fetcher = None
    web_cache = WebCache(fetcher)
    webcache_path = os.path.join(os.getcwd(), WEBCACHE_PATH)
    if pytestconfig.option.webcache_load:
        web_cache.load_from_dir(webcache_path)

    yield web_cache  # each test runs here

    # --- teardown ---
    web_cache.report_stats()
    if pytestconfig.option.webcache_save:
        web_cache.save_to_dir(webcache_path)


@pytest.fixture(autouse=True)
def cache_url_data_for_plugin(monkeypatch, web_cache):
    """Cache data returned from URLs fetched over the network.

    Wrap URL fetcher cache around code so we don't hit the real PDB with
    repeated requests from this test.

    Note that his fixture can't be session scoped like the cache itself since
    the monkeypatch dependence fixture is function scoped.
    """
    # --- setup ---
    monkeypatch.setattr(plugin.pdb._fetcher, 'get_data', web_cache.access)

    yield  # each test runs here

    # --- teardown ---
    pass  # monkeypatch undoes all recorded changes


@pytest.fixture(autouse=True)
def chdir_back_to_pwd(monkeypatch):
    """Make sure we cd back to PWD after test completes (failure or not)."""
    # --- setup ---
    monkeypatch.chdir('.')

    yield  # each test runs here

    # --- teardown ---
    pass  # monkeypatch undoes all recorded changes


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
            return json.dumps(self._data)

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

    # When unset (first time use of plugin), set preference to WARNING.
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


def test_run_main_pdbid_file_exists():
    """Load and initialize the module and attempt calling main().
    This test covers the special case where a pdbid was requests but we happen
    to also find the corresponding .cif file locally in the working directory.
    """
    pdbid = '3mxw'
    os.chdir('tests/data')
    assert os.access(pdbid + '.cif', os.R_OK)
    plugin.main([pdbid])


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


def test_chimera_molecule():
    """Tests PDB entry with > 1 molecule_name."""
    plugin.PDB_Analysis_Molecules('1f0d')
    assert plugin.count_chains() == 1


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


# ----- keep this test last -----
def test_print_web_cache_stats(web_cache, capsys):
    """Not a test: Just print web_cache stats outside of stdout capture."""
    with capsys.disabled():
        web_cache.report_stats()
    assert True
