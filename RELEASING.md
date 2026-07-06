# Releasing

Releases are published to PyPI automatically when a GitHub Release is created —
there is no manual upload step and no stored API token.

## Cut a release

1. Make sure `main` is green and contains everything for the release.
2. Create a GitHub Release with a semantic-version tag (the leading `v` is
   stripped by `setuptools_scm`, so `v0.2.0` becomes version `0.2.0`):

   ```bash
   gh release create v0.2.0 --generate-notes
   ```

   `--generate-notes` fills the body from the merged pull requests, grouped by
   label according to [`.github/release.yml`](.github/release.yml). Review and
   edit the notes, then publish the release.
3. Publishing the release triggers
   [`.github/workflows/publish.yml`](.github/workflows/publish.yml), which builds
   the sdist + wheel and uploads them to PyPI via
   [Trusted Publishing](https://docs.pypi.org/trusted-publishers/) (OIDC).

## Versioning

- The version is derived from the git tag by `setuptools_scm`; never hardcode a
  version in `pyproject.toml`.
- Between releases the version is a development version (e.g.
  `0.2.0.dev5+g<sha>`), a PEP 440 *local* version that PyPI rejects — so only a
  tagged release can be published (a built-in safety net against accidental
  uploads).

## Notes

- Label pull requests (`enhancement`, `bug`, `documentation`, `dependencies`,
  ...) so the generated notes are grouped; unlabelled PRs still appear under
  **Other Changes**.
- The PyPI trusted publisher is configured for the project `thermoscreening`
  (owner `MolarVerse`, repository `ThermoScreening`, workflow `publish.yml`, no
  environment). If a release's publish job fails with `invalid-publisher`, these
  four values on PyPI must match the workflow.
- To publish **conda-forge** as well (the channel that also pulls the DFTB+ /
  xtb / tblite backends), `PQAnalysis` first needs to be available on
  conda-forge; then submit a recipe via `staged-recipes`.
