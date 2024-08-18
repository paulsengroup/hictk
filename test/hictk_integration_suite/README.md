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
venv/bin/hictk_integration_suite /tmp/hictk/build/src/hictk/hictk config.toml --data-dir=../data

{'title': 'hictk-metadata-cli', 'args': ['/tmp/hictk/build/src/hictk/hictk', 'metadata'], 'elapsed-time': datetime.timedelta(microseconds=14914), 'exit-code': 106, 'notes': [], 'status': 'PASS'}
{'title': 'hictk-metadata-cli', 'args': ['/tmp/hictk/build/src/hictk/hictk', 'metadata', '--help'], 'elapsed-time': datetime.timedelta(microseconds=10365), 'exit-code': 0, 'notes': [], 'status': 'PASS'}
{'title': 'hictk-metadata-cli', 'args': ['/tmp/hictk/build/src/hictk/hictk', 'metadata', 'not-a-file'], 'elapsed-time': datetime.timedelta(microseconds=9756), 'exit-code': 105, 'notes': [], 'status': 'PASS'}
{'title': 'hictk-metadata-cli', 'args': ['/tmp/hictk/build/src/hictk/hictk', 'metadata', '--foobar'], 'elapsed-time': datetime.timedelta(microseconds=9833), 'exit-code': 106, 'notes': [], 'status': 'PASS'}
...

# PASS: 35
# FAIL: 0
```
