# Data

| | |
|---|---|
| `sample/` | Complete output for four models (`ellipse`, `chair_remesh`, `vase`, `kitten`), tracked in this repository so the method's output can be inspected without a download. |
| everything else | Published as release assets — see below. |

The input meshes for `ellipse` and `chair_remesh` are in `../test_data/`, so
those two can be re-run end to end straight from a clone.

## Full dataset

Results for all 125 models, and the input meshes that produced them, are
attached to the [latest release](https://github.com/nicopietroni/LoopWeaver/releases/latest):

| Asset | Size | Contents |
|---|---|---|
| `flowloops-inputs.tar.gz` | 50 MB | the 125 input triangle meshes |
| `flowloops-results.tar.gz` | 598 MB | curve networks, surfaced meshes, cross fields, timings, screenshots |

```bash
curl -LO https://github.com/nicopietroni/LoopWeaver/releases/latest/download/flowloops-inputs.tar.gz
curl -LO https://github.com/nicopietroni/LoopWeaver/releases/latest/download/flowloops-results.tar.gz
tar xzf flowloops-inputs.tar.gz
tar xzf flowloops-results.tar.gz
```

Each archive contains its own `README.md` describing the per-model file layout.

## Which version produced these results

The results in `sample/`, `output/` and the release archives were produced by
an earlier version of LoopWeaver: the code at the `v1.0` and `pre-nico-update`
tags. Re-running the current version will not reproduce them exactly. Its
defaults differ (sampling density 6000 instead of 10000, maximum side angle
120° instead of 180°), the batch scripts add `-error_setup 1 -subsampl 2`,
and each archived run used its own reconstruction error (0.02 for 101
models, 0.03 for 16, 0.04 for 8; see each model's `_config.txt`). With
symmetry, the archived `.path` files may also hold fewer loops than the
current ones, which cover the whole mirrored mesh.
