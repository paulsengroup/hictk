<!--
Copyright (C) 2024 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# README

Installing the pipeline:

```bash
python3 -m venv venv
venv/bin/pip install . -v
```

Running the tests:

```console
user@dev:/tmp$ venv/bin/hictk_integration_suite /tmp/hictk/build/src/hictk/hictk config.toml --data-dir=../data
INFO:root:staging test files for balance tests took 18.58ms
INFO:root:planning for balance tests took 86.97ms
INFO:root:{'id': '91501f4f7cb8da7eb1f26adcdc98681af27f88770c959531325f244fe0fbf8f0', 'title': 'hictk-balance-ice-cli', 'args': ['balance', 'ice'], 'hictk-runtime': '0:00:00.004082', 'validation-runtime': '0:00:00.000006', 'elapsed-time': '0:00:00.004088', 'exit-code': 106, 'expect-failure': True, 'errors': [], 'status': 'PASS'}
INFO:root:{'id': '78fdc5f107b74d28b25e759984f1fdd74a2f0c070f20ca5606a736200f2d1c9f', 'title': 'hictk-balance-ice-cli', 'args': ['balance', '/tmp/hictk-integration-test-t3g8x7bu/staged_files/ENCFF993FGR.2500000.cool'], 'hictk-runtime': '0:00:00.003897', 'validation-runtime': '0:00:00.000005', 'elapsed-time': '0:00:00.003902', 'exit-code': 106, 'expect-failure': True, 'errors': [], 'status': 'PASS'}
INFO:root:{'id': '732a3fb663deb58b41dd8c18dbd559a67db6ddb29ee6eda51e667a7e6d81b606', 'title': 'hictk-balance-ice-cli', 'args': ['balance', 'ice', '/tmp/hictk-integration-test-t3g8x7bu/staged_files/ENCFF993FGR.2500000.cool', '--mode=invalid'], 'hictk-runtime': '0:00:00.003845', 'validation-runtime': '0:00:00.000003', 'elapsed-time': '0:00:00.003847', 'exit-code': 105, 'expect-failure': True, 'errors': [], 'status': 'PASS'}
...
INFO:root:running tests for zoomify took 15.45s
INFO:root:running 546 tests took 247.72s

# PASS: 546
# SKIP: 0
# FAIL: 0
```
