Real Gaussian 16 log files, used to test `read_cclib`/`cclib_thermo` against
genuine (non-mocked) Gaussian output.

- `water_neutral_opt_freq.out` — a clean neutral water optimization + frequency
  job. Parses successfully; used as the real-data success-path test.
- `mp2_avdz_freq_tight.log` — an MP2/aug-cc-pVDZ frequency job whose "Leave
  Link" timing line (`MaxMem=... cpu: ... elap: ...`) crashes cclib 1.8.1's
  Gaussian parser (confirmed: still the latest cclib release at the time).
  Root cause is an unhandled `ValueError` inside cclib's own parser; what
  actually reaches `read_cclib` depends on the surrounding logging setup (a
  bare `cclib.io.ccread` call raises that `ValueError` directly, but cclib's
  own `logger.error()` call right before it can itself raise first under some
  logging configurations -- either way, `read_cclib` wraps whatever comes out
  into a clean `TSValueError`). Used as the regression test for that
  wrapping, so no cclib parsing exception escapes uncontrolled.

Source: [cclib/cclib-data](https://github.com/cclib/cclib-data), the cclib
project's own public regression-test data, `Gaussian/Gaussian16/`.
