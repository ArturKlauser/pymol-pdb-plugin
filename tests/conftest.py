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
