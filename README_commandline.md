# LoopWeaver — command line

```
loop_weaver mesh.obj [mesh.field] [options]
```

Without `-batch` the GUI opens with the mesh loaded and the options applied.
With `-batch` the whole pipeline runs (symmetry cut, cleaning, cross field,
decomposition, relaxation, final reconstruction), every output is saved next
to the mesh, and the program quits.

Running `loop_weaver` without arguments prints the list of options.

## Parameters: settings file and command line

The parameters are set in three layers, each overriding the previous one:

1. the compiled defaults;
2. the settings file `loopweaver_settings.txt` next to the executable. If it
   is missing it is written with the compiled defaults; delete it to write it
   again (for example after an update that changed the defaults). Every
   parameter of the GUI is in this file, one per line, `name value  # comment`;
3. the command line options below.

`-settings <file>` reads another settings file instead of
`loopweaver_settings.txt`. The `<name>_settings.txt` saved with every result
has the same format, so a run can be repeated with
`-settings mesh_settings.txt`.

The `-x` / `-no_x` pairs switch an option on or off regardless of the
settings file.

## Options

### Input and run

| Option | Description |
|---|---|
| `mesh.obj` | input triangle mesh (non-manifold elements are split on load) |
| `mesh.field` | optional cross field of the mesh (given after the mesh); without it the field is computed from the curvature |
| `-batch` | process, save everything and quit |
| `-settings <file>` | settings file to use instead of `loopweaver_settings.txt` |

### Decomposition

| Option | Default | Description |
|---|---|---|
| `-symm` / `-no_symm` | off | cut the mesh with the symmetry plane x=0, decompose one half and mirror the result |
| `-angle <deg>` | 50 | sharp features angle (dihedral); 180: no sharp features |
| `-split` / `-no_split` | off | split the paths at the crossings when removing them at the end |
| `-complete_features` / `-no_complete_features` | off | complete all the sharp features with paths |
| `-dynamic` / `-no_dynamic` | on | dynamic update: retrace the candidates after every insertion |
| `-relax` / `-no_relax` | on | relaxation at the end: moves the paths to spread them (graph distances from the other paths, the borders and the sharp features), then tries to remove paths and relaxes again |

### Patch conditions

| Option | Default | Description |
|---|---|---|
| `-error <ratio>` | 0.02 | loop reconstruction condition, max distance of the reconstructed patch (× bounding box diagonal); also switches the condition on |
| `-no_error` | | no loop reconstruction condition |
| `-normal_error <deg>` | 30 | max normal error of the reconstruction (on a fraction of the faces, `normal_angle_error_percent`, 0.5); also switches it on |
| `-no_normal` | | reconstruction checked on the distance only |
| `-error_setup <0..3>` | | preset of the reconstruction errors (as the old version): 0 = 0.02 / 30° / 0.2, 1 = 0.03 / 40° / 0.3, 2 = 0.04 / 50° / 0.5, 3 = 0.05 / 60° / 0.6 (distance ratio / normal angle / percentile); a negative value switches the condition off |
| `-side_angle <deg>` | 120 | max turning of a patch side (sum of the 3D angles along the side); also switches the condition on |
| `-no_side` | | no side curvature condition |
| `-side_smooth <n>` | 10 | Taubin steps on the side polyline before measuring its angle (removes the zig-zag along the mesh edges, keeps the real curvature); 0: none |
| `-bezier <percent>` | 2 | max error of the quadratic Bézier fitted on each side (percent of the bounding box diagonal); also switches the condition on |
| `-no_bezier` | | no Bézier error condition |
| `-gauss <deg>` | off (90) | developability: max total absolute Gaussian curvature of a patch (sum of the angle defects of its interior vertices); also switches the condition on |
| `-no_gauss` | | no Gaussian curvature condition |
| `-gauss_smooth <n>` | 5 | diffusion steps of the signed angle defect before taking the absolute value (removes the tessellation noise) |

### Surfaces

| Option | Default | Description |
|---|---|---|
| `-thinplate` / `-no_thinplate` | off | thin-plate only surfaces, both in the reconstruction condition and in the final reconstruction |
| `-subsampl <n>` | 1 | subsampling of the paths in the reconstruction condition |

Unknown options are reported and ignored.

## Output files

With `-batch` (or *Save All* in the GUI), next to the mesh `name.obj`:

| File | Content |
|---|---|
| `name.path` | the loops, as sequences of vertex indices of `name_remeshed.obj`. Each path starts with a header line `P <i>`; the paths lying on the sharp features are marked `P <i> S` |
| `name_remeshed.obj` | the processed mesh the paths refer to (cleaned, refined; with symmetry the mirrored whole mesh) |
| `name_edge.obj` | the loops as an edge mesh |
| `name_solved.obj` | the final reconstructed surface |
| `name.field` | the cross field used |
| `name_timings.txt` | field, tracing and surfacing times, and the total length of the loops |
| `name_settings.txt` | all the parameters of the run (settings file format, can be passed back with `-settings`) |

The exchange files of the surfacing solver (`.curve`, `.cycle`,
`.normalcurve`, `.feat`, `.curve_types`) are deleted at the end.

## Examples

```bash
# all the defaults (settings file)
./loop_weaver model.obj -batch

# symmetric model, sharp features at 45 degrees, paths subsampled
./loop_weaver model.obj -batch -symm -angle 45 -subsampl 2

# no sharp features, no dynamic update, thin-plate surfaces
./loop_weaver model.obj -batch -angle 180 -no_dynamic -thinplate

# looser reconstruction, developable patches
./loop_weaver model.obj -batch -error_setup 2 -gauss 60

# repeat a previous run
./loop_weaver model.obj -batch -settings model_settings.txt
```

## Batch scripts

The four scripts `run_weaver_{features,no_features}_{symm,no_symm}.sh` run
`loop_weaver -batch` on all the `.obj` of a folder:

```bash
./run_weaver_features_symm.sh <folder> [extra loop_weaver options]
```

| Script | Options |
|---|---|
| `run_weaver_features_symm.sh` | `-symm -angle 45 -subsampl 2` |
| `run_weaver_features_no_symm.sh` | `-no_symm -angle 45 -subsampl 2` |
| `run_weaver_no_features_symm.sh` | `-symm -angle 180 -subsampl 2` |
| `run_weaver_no_features_no_symm.sh` | `-no_symm -angle 180 -subsampl 2` |

The extra options after the folder are passed to every run (for example
`-no_dynamic -thinplate`). Each mesh has a timeout of 1800 s; the processed
meshes are moved to `<folder>/done`, and every run (options, success,
error, timeout, time) is logged in `timeouts.txt`.
