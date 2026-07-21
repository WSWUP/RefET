import os

import pytest


def pytest_collection_modifyitems(config, items):
    """Skip tests marked 'network' unless RUN_NETWORK_TESTS is set.

    The network tests hit the live AgriMet service; everything else runs
    offline against recorded data.
    """
    if os.environ.get('RUN_NETWORK_TESTS'):
        return
    skip = pytest.mark.skip(
        reason='network test; set RUN_NETWORK_TESTS=1 to run')
    for item in items:
        if 'network' in item.keywords:
            item.add_marker(skip)
