import pymol
import pytest


@pytest.hookimpl(hookwrapper=True)
def pytest_terminal_summary(terminalreporter, exitstatus, config):
    yield
    if exitstatus:
        print('pytest encountered errors')
    print('PyMOL shutting down with exitstatus %d' % exitstatus)
    # Under python2 we don't get non-0 exit status propagated out of a failing
    # pytest run. The cuplrit appears to lie in something that happens in the
    # pymol.finish_lauching() starup code, which is swallowing pytest's exit
    # status and always returns 0 instead.  This is a problem for the CI
    # workflow as it doesn't detect that tests have failed and thus doesn't
    # notify us. Explicitly calling PyMOL's quit() function works around this
    # issue.
    pymol.cmd.quit(exitstatus)


def add_bool_option(parser, name, default=False, help='', **kwargs):
    """Adds a boolean option with --name and --no-name variants."""
    positive_option = '--' + name
    negative_option = '--no-' + name
    dest = name.replace('-', '_')  # make it a proper python name
    # Remove keys from kwargs which we set ourselves.
    for key in ('dest', 'action', 'nargs'):
        kwargs.pop(key, None)
    # Add --name option.
    parser.addoption(
        positive_option,
        dest=dest,
        action='store_true',
        # nargs=0,
        default=default,
        help='%s (default: %s)' % (help, default),
        **kwargs)
    # Add --no-name option.
    parser.addoption(
        negative_option,
        dest=dest,
        action='store_false',
        # nargs=0,
        default=not default,
        help='Turn off ' + positive_option,
        **kwargs)


def pytest_addoption(parser):
    # Webcache options:
    #
    # This test suite uses a cache to speed up web accesses by the plugin and
    # PyMOL. By default it is configured to run entirely isolated from the
    # network, which also speeds up testing significantly. The options have the
    # following meaning:
    #   --[no-]webcache: (default: on)
    #     Turn web cache on/off. If off, the original access paths are used (no
    #     patching happens) and all other webcache* options are meaningless.
    #   --[no-]webcache-load: (default: on)
    #     Web cache data is loaded from the cache directory on session start.
    #     Used to avoid network accesses for data that has been recorded and
    #     saved in the cache directory.
    #   --[no-]webcache-fetch: (default: off)
    #     Cache misses are fetched over the network via the original access
    #     paths. If off, any cache miss results in an exception.
    #   --[no-]webcache-save: (default: off)
    #     Web cache data is saved to the cache directory on session end.
    #
    # Common operations:
    #   Running tests:
    #     Use the default option settings.
    #   Changing tests or writing new tests:
    #     Chances are that you are affecting the accesses necessary to run the
    #     tests. While debugging, to fetch any missing data, turn on
    #       --webcache-fetch
    #   Committing changed tests:
    #     To repopulate the cache directory perform a passing run with
    #       --no-webcache-load --webcache-fetch --webcache-save
    #     Then git add/rm the respective changes in the cache directory.

    add_bool_option(parser,
                    'webcache',
                    default=True,
                    help='Use a web cache to filter network traffic.')
    add_bool_option(parser,
                    'webcache-load',
                    default=True,
                    help='Load web cache from directory before runs start.')
    add_bool_option(parser,
                    'webcache-fetch',
                    default=False,
                    help='Fetch items missing in web cache from web.')
    add_bool_option(parser,
                    'webcache-save',
                    default=False,
                    help='Save web cache to directory after runs complete.')
