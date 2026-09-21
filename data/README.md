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
