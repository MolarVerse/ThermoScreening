# conda-forge recipes

Draft [conda-forge](https://conda-forge.org/) recipes for distributing
ThermoScreening (and its dependency PQAnalysis) through the `conda-forge`
channel. They are kept here for reference and maintenance; the recipes that
conda-forge actually builds live in per-package *feedstock* repositories created
from [`conda-forge/staged-recipes`](https://github.com/conda-forge/staged-recipes).

## Why two recipes

`ThermoScreening` depends on `PQAnalysis`, which is **not yet on conda-forge**.
conda-forge packages may only depend on other conda-forge packages, so
`PQAnalysis` has to land first. Every other dependency (`numpy`, `scipy`,
`pymatgen-core`, `beartype`, `ase`, `rdkit`, and PQAnalysis's own
`multimethod`/`lark`/`tqdm`/`decorator`/`argcomplete`/`rich-argparse`) is
already available on conda-forge.

| Recipe | noarch? | Notes |
|--------|---------|-------|
| `pqanalysis/` | no | Ships a compiled Cython extension, so it builds per platform (needs a C compiler). |
| `thermoscreening/` | yes | Pure Python. |

## Submission order

1. **PQAnalysis first.** Fork `conda-forge/staged-recipes`, copy
   `pqanalysis/` into its `recipes/` directory, and open a PR. Once it is
   merged, conda-forge's bot creates `PQAnalysis-feedstock` and publishes the
   package (usually within an hour).
2. **ThermoScreening second.** After `pqanalysis` is available on the
   `conda-forge` channel, submit `thermoscreening/` the same way. Its
   `pqanalysis >=1.3.0` run requirement will then resolve.

staged-recipes can build sibling recipes in dependency order within a single
PR, so submitting both at once can work — but the two-step order above is the
simpler, lower-risk path and lets `pqanalysis` publish before `thermoscreening`
is reviewed.

## Before submitting

- **Maintainer(s):** `extra.recipe-maintainers` lists `galjos`. Add any other
  GitHub usernames who should co-maintain the feedstocks.
- **Versions & hashes** are pinned to the current PyPI releases
  (PQAnalysis 1.3.0, ThermoScreening 0.1.0). To refresh for a new release,
  bump `version` and replace `sha256` with the sdist hash:

  ```bash
  # prints the sha256 of the PyPI source tarball
  curl -sL https://pypi.org/pypi/ThermoScreening/json \
    | python -c "import json,sys; d=json.load(sys.stdin); \
        print(next(u['digests']['sha256'] for u in d['urls'] if u['packagetype']=='sdist'))"
  ```

  After the feedstocks exist, conda-forge's `regro-cf-autotick-bot` opens
  version-bump PRs automatically, so this is mainly needed for the initial
  submission.

## Local check (optional)

If you have `conda-build` installed you can lint/build a recipe before
submitting:

```bash
conda build conda-recipes/pqanalysis
conda build conda-recipes/thermoscreening -c conda-forge
```

See the conda-forge [contributing guide](https://conda-forge.org/docs/maintainer/adding_pkgs/)
for the full staged-recipes workflow.
