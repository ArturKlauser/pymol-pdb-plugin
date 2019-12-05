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
    add_bool_option(parser,
                    'webcache-save',
                    default=False,
                    help='Save web cache to directory after runs complete.')
    add_bool_option(parser,
                    'webcache-load',
                    default=True,
                    help='Load web cache from directory before runs start.')
    # By default we isolate the test from network fetch dependencies.
    add_bool_option(parser,
                    'webcache-fetch',
                    default=False,
                    help='Fetch items missing in web cache from web.')
