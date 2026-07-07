# conda-forge recipe

Draft [conda-forge](https://conda-forge.org/) recipe for distributing
ThermoScreening through the `conda-forge` channel. It is kept here for reference
and maintenance; the recipe conda-forge actually builds lives in the
`thermoscreening` *feedstock* created from
[`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes).

| Recipe | noarch? | Notes |
|--------|---------|-------|
| `thermoscreening/` | yes | Pure Python. |

## Dependency: PQAnalysis

`ThermoScreening` depends on `PQAnalysis`, which must be on conda-forge first
(conda-forge packages may only depend on other conda-forge packages). Its recipe
lives in the [PQAnalysis repository](https://github.com/MolarVerse/PQAnalysis)
(`conda-recipes/pqanalysis/`). Every other dependency (`numpy`, `scipy`,
`pymatgen-core`, `beartype`, `ase`, `rdkit`) is already on conda-forge.

Both recipes were submitted together to `staged-recipes`, which builds sibling
recipes in dependency order (`pqanalysis` first, then `thermoscreening`).

## Before submitting

- **Maintainer(s):** `extra.recipe-maintainers` lists `galjos`. Add any other
  GitHub usernames who should co-maintain the feedstock.
- **Versions & hashes** are pinned to the current PyPI release
  (ThermoScreening 0.1.0). To refresh for a new release, bump `version` and
  replace `sha256` with the sdist hash:

  ```bash
  # prints the sha256 of the PyPI source tarball
  curl -sL https://pypi.org/pypi/ThermoScreening/json \
    | python -c "import json,sys; d=json.load(sys.stdin); \
        print(next(u['digests']['sha256'] for u in d['urls'] if u['packagetype']=='sdist'))"
  ```

  After the feedstock exists, conda-forge's `regro-cf-autotick-bot` opens
  version-bump PRs automatically, so this is mainly needed for the initial
  submission.

## Local check (optional)

If you have `conda-build` installed you can lint/build the recipe before
submitting:

```bash
conda smithy recipe-lint conda-recipes/thermoscreening
conda build conda-recipes/thermoscreening -c conda-forge
```

See the conda-forge [contributing guide](https://conda-forge.org/docs/maintainer/adding_pkgs/)
for the full staged-recipes workflow.
