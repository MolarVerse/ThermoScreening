Real Gaussian 16 log files, used to test `read_cclib`/`cclib_thermo` against
genuine (non-mocked) Gaussian output.

- `water_neutral_opt_freq.out` — a clean neutral water optimization + frequency
  job. Parses successfully; used as the real-data success-path test.
- `mp2_avdz_freq_tight.log` — an MP2/aug-cc-pVDZ frequency job whose "Leave
  Link" timing line (`MaxMem=... cpu: ... elap: ...`) crashes cclib 1.8.1's
  Gaussian parser with an unhandled `ValueError` (confirmed: still the latest
  cclib release at the time). Used as the regression test for wrapping any
  cclib parsing exception into a clean `TSValueError` instead of letting it
  escape uncontrolled.

Source: [cclib/cclib-data](https://github.com/cclib/cclib-data), the cclib
project's own public regression-test data, `Gaussian/Gaussian16/`.
