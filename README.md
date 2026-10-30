# A phase gap between transcriptional peaks of clock neurons encodes seasonal time

Analysis code accompanying the study, covering bootstrap statistics, kernel-based activity traces, lipid-droplet quantification, RNAscope quantification, and live-imaging cell tracking and MS2 signal processing.

The MATLAB figure analyses include their input spreadsheets. The image-analysis workflows require user-supplied raw images and manual annotations. Three initial model files are included for cell tracking; per-recording fine-tuning produces the weights used for subsequent tracking.

## Repository contents

| Folder | Purpose | Entry point |
| --- | --- | --- |
| [Bootstrap analysis](Bootstrap%20analysis/) | Bootstrap estimates of peak times, phase gaps, peak shifts, rhythmicity and amplitude | One MATLAB script per figure folder |
| [KDE analysis](KDE%20analysis/) | Binary activity KDE, mean-activity dot plots, continuous recordings and light-pulse activity curves | One MATLAB script per figure folder |
| [Lipid droplet analysis](Lipid%20droplet%20analysis/) | Guided central-brain ROI preparation, droplet detection, cross-Z linking, QC and summaries | [Lipid_droplet_workflow.ipynb](Lipid%20droplet%20analysis/Lipid_droplet_workflow.ipynb) |
| [RNAscope analysis](RNAscope%20analysis/) | LD4 single-channel and LP dual-channel spot quantification | [RNAscope_quantification_usage.ipynb](RNAscope%20analysis/RNAscope_quantification_usage.ipynb) |
| [Cell_tracking](Cell_tracking/) | First-volume segmentation, manual correction, U-Net fine-tuning and 3D cell tracking | [Cell_tracking.ipynb](Cell_tracking/Cell_tracking.ipynb) |
| [Signal_processing](Signal_processing/) | MS2 signal extraction, manual review, baseline correction and numerical/figure export | [Signal_processing.m](Signal_processing/Signal_processing.m) |

Each module keeps its shared functions in its own `script/` directory, except `Cell_tracking`, which provides replacement files for an upstream package. Keep these directory structures when moving the repository. Folder and script names follow the current manuscript figure numbering.

In the M/E analyses, **M** denotes sLNv cells and **E** denotes the 5th LN. Fig1H analyzes lLNv cells. Time units are specified in the input headers and configuration; recording timestamps are generally stored in hours, while kernel windows are configured in minutes.

## Software requirements

### MATLAB

The MATLAB workflows were prepared and checked with **MATLAB R2023a**.

| Module | Additional products |
| --- | --- |
| Bootstrap analysis | Statistics and Machine Learning Toolbox |
| KDE analysis | Statistics and Machine Learning Toolbox for the complete module, including Fig1H statistics |
| Signal_processing | Image Processing Toolbox, Statistics and Machine Learning Toolbox, and Curve Fitting Toolbox |
| Signal_processing, optional PPT export | MATLAB Report Generator |

Cell identification and signal proofreading require the MATLAB desktop and interactive figure windows.

### Python image quantification

The lipid-droplet and RNAscope modules use **Python 3.11.9** as their reference version. Each contains a `requirements.txt` with pinned package versions. They can share a Python 3.11 environment. For example, from the repository root in an Anaconda PowerShell Prompt:

```powershell
conda create -n phasegap_imaging python=3.11.9
conda activate phasegap_imaging
python -m pip install -r "Lipid droplet analysis/requirements.txt" -r "RNAscope analysis/requirements.txt"
python -m pip install notebook
python -m ipykernel install --user --name phasegap_imaging --display-name "Phase-gap image analysis"
```

Open the desired notebook and select that kernel:

```powershell
python -m notebook "Lipid droplet analysis/Lipid_droplet_workflow.ipynb"
```

For RNAscope, open `RNAscope analysis/RNAscope_quantification_usage.ipynb` instead. Both notebooks locate their module from the working directory or its parents. Set `PROJECT_ROOT` explicitly if the kernel starts elsewhere.

**Cell tracking uses a separate historical Win10 environment**, described below. Its TensorFlow/Keras environment is separate from the Python 3.11 quantification environment.

## Bootstrap analysis

### Run a figure

Open the figure's saved `.m` file and use **Run/F5**, or run it by its full path:

```matlab
repoDir = 'C:/path/to/this/repository';
run(fullfile(repoDir, 'Bootstrap analysis', 'Fig2F', 'Fig2F.m'));
```

The entry script locates its input files and adds the module's shared `script/` folder automatically. Running the complete saved file supplies the script location needed by `mfilename('fullpath')`.

### Figure entry points

| Entry script | Analysis |
| --- | --- |
| [Fig2C](Bootstrap%20analysis/Fig2C/Fig2C.m) | Fly-paired KDE and continuous-intensity analyses |
| [Fig2F](Bootstrap%20analysis/Fig2F/Fig2F.m) | Paired M/E peak-time and phase-gap bootstrap |
| [Fig2H](Bootstrap%20analysis/Fig2H/Fig2H.m) | Paired M/E peak-time and phase-gap bootstrap |
| [Fig5B](Bootstrap%20analysis/Fig5B/Fig5B.m) | Paired M/E peak-time and phase-gap bootstrap |
| [Fig5C](Bootstrap%20analysis/Fig5C/Fig5C.m) | Paired M/E peak-time and phase-gap bootstrap |
| [Fig5DE](Bootstrap%20analysis/Fig5DE/Fig5DE.m) | M- and E-cell peak shifts relative to the control |
| [Fig5F](Bootstrap%20analysis/Fig5F/Fig5F.m) | Paired M/E peak-time and phase-gap bootstrap |
| [Fig5H](Bootstrap%20analysis/Fig5H/Fig5H.m) | Both dbtS and dbtL phase-gap contrasts against their respective controls |
| [figS5C](Bootstrap%20analysis/figS5C/figS5C.m) | Paired M/E peak-time and phase-gap bootstrap |
| [figS5E](Bootstrap%20analysis/figS5E/figS5E.m) | Periodic-spline rhythmicity, amplitude and equivalence analyses |
| [figS6C](Bootstrap%20analysis/figS6C/figS6C.m) | Whole-cell and nuclear measurements, run separately |
| [figS7C](Bootstrap%20analysis/figS7C/figS7C.m) | Paired M/E peak-time and phase-gap bootstrap |
| [figS8A](Bootstrap%20analysis/figS8A/figS8A.m) | Peak shifts for three experimental groups relative to the control |
| [figS10](Bootstrap%20analysis/figS10/figS10.m) | Paired analyses across the configured photoperiod datasets |
| [figS11B](Bootstrap%20analysis/figS11B/figS11B.m) | Paired M/E peak-time and phase-gap bootstrap |
| [figS13CFI](Bootstrap%20analysis/figS13CFI/figS13CFI.m) | Separate analyses for supplementary panels C, F and I |

### Inputs and settings

Most figure folders contain `data.xlsx`. The input blocks are selected by sheet index and, where needed, explicit cell ranges in the entry script. Preserve tab order, column positions and paired M/E rows. For these index-based readers, changing a tab's display name does not change the selected data.

**Fig2C uses two workbooks**:

- `data_binary.xlsx`: binary ON/OFF activity used for KDE.
- `data_intensity.xlsx`: continuous signal values used for the intensity branch.

Each Fig2C tab represents one recording and uses a `YYMMDD FLYn` name. Corresponding tabs are matched by this date/animal key across the two workbooks. Each cell occupies adjacent `<cell>_Time_h` and `<cell>_Value` columns. Nominal start CT is inferred by rounding the first valid timestamp, without wrapping CT24 to zero. Binary-only recordings remain available to the KDE branch; intensity analysis uses matched recordings. The current batch split is determined by the binary tab order and `batchSheetCounts = [38, 15]`.

The figure scripts specify the bootstrap counts, random seeds, effect-size margins, significance level and display settings. Keep the supplied analysis settings when reproducing the supplied analyses. Figure dimensions and plotting options can be adjusted independently where indicated in the scripts.

Contrast signs are explicit in the entry scripts and outputs. The standard single-group phase-gap scripts report **M minus E**. Fig5H reports **(E minus M in the condition) minus (E minus M in the control)**. Peak-shift comparisons report the experimental peak minus the control peak for each cell type.

### Outputs

Outputs are written to the figure's `results/` folder. Depending on the analysis, they include XLSX summary and diagnostic tables, MAT results with effective settings, and SVG plots. Fig2C saves its KDE and intensity branches under separate output stems. Change the configured output directory to retain separate repeated runs. Figure visibility is controlled by `showFigures` in the relevant configuration block.

## KDE analysis

Run these scripts in the same way as the bootstrap entry points:

```matlab
repoDir = 'C:/path/to/this/repository';
run(fullfile(repoDir, 'KDE analysis', 'Fig2B', 'Fig2B.m'));
```

| Entry script | Analysis | Kernel window |
| --- | --- | --- |
| [Fig1F](KDE%20analysis/Fig1F/Fig1F.m) | Aligned traces at ZT4, ZT8, ZT12 and ZT24 | 1 min |
| [Fig1H](KDE%20analysis/Fig1H/Fig1H.m) | Method7/7b KDE-mean jitter dot plots, cell- and fly-level summaries, and fly-level Kruskal–Wallis comparison | 20 min |
| [Fig2B](KDE%20analysis/Fig2B/Fig2B.m) | Individual M/E traces at the configured CT groups | 20 min |
| [Fig2D](KDE%20analysis/Fig2D/Fig2D.m) | Continuous recordings, including records longer than 24 h | 180 min |
| [Fig3CFI](KDE%20analysis/Fig3CFI/Fig3CFI.m) | M/E light-pulse curves and pre/post comparisons at ZT15, ZT18 and ZT21 | 30 min |

Fig1F, Fig1H, Fig2B and Fig2D read `data_binary.xlsx`. Their tabs contain adjacent `<cell>_Time_h` / `<cell>_Value` pairs, with a separate native time axis for each trace. Preserve tab and column order. In particular, Fig2B tab names retain their CT labels, and the first eight Fig2D tabs retain their legacy shared-time interpretation. Recording tabs use date/animal labels where those identifiers were available.

`script/compute_KDE.m` applies the original Gaussian kernel to binary activity, with standard deviation equal to the configured window divided by six. NaN signal values are treated as zero by this function. In the distributed Fig2B table, blank signal values at recorded time points have been explicitly filled with zero.

Fig2D exports full-length traces. Its summary uses the original absolute-time selection `t < 24 h`, pools cells, and retains the original zero-fill and smoothing rules for the summary. Recording times above 24 h remain in the full trace outputs.

**Fig3CFI uses a different estimator.** Its `data.xlsx` contains fractional activity matrices, with time in the first column and one animal per following column. Tabs retain names such as `EcellZT15`, `EcellZT15_control`, `McellZT15` and `McellZT15_control`, repeated for ZT18 and ZT21. Keep these names. `compute_NW_activity.m` performs locally normalized Nadaraya–Watson Gaussian regression and preserves fractional values. Pre/post statistical comparisons use the raw observations, with three pre-LP and two post-LP time points. The script exports the combined M/E figures and the individual C/F/I panels.

Results are written to each figure's `results/` folder: numerical tables, MAT results, parameter JSON files and figures. `cfg.showFigures` controls window visibility, and `cfg.exportPNG` controls PNG export. Input QA reports preserve information about missing values and repeated or nonmonotonic timestamps.

## Lipid-droplet quantification

Open [Lipid_droplet_workflow.ipynb](Lipid%20droplet%20analysis/Lipid_droplet_workflow.ipynb) in the Python 3.11 environment.

### Configure the inputs

Set `DATA_ROOT` in the notebook to the parent directory containing the IMS files or batches. Edit [config/workflow.json](Lipid%20droplet%20analysis/config/workflow.json) to select `input_paths`, output location, group labels, group order and statistical comparisons. Edit [config/analysis.yaml](Lipid%20droplet%20analysis/config/analysis.yaml) for ROI and detection parameters when adapting the workflow to a new acquisition.

Paths in the JSON configuration are resolved relative to that configuration file; entries in `input_paths` are resolved under `data_root`. The notebook's `DATA_ROOT` value overrides the configured data root. The referenced `example_data/` location is a placeholder for user-supplied images.

### Run the notebook in stages

1. Inspect the guide-export plan. Enable `RUN_EXPORT` to export the native-intensity, level-0 maximum projections. This creates `central_brain_roi/` beside the input IMS files and writes `guide_manifest.json`.
2. Open each exported projection in Fiji, draw one closed central-brain area ROI, and save it with the projection's basename and a `.roi` extension. Keep the projection's original geometry while annotating.
3. Inspect the quantification plan, then enable `RUN_QUANTIFICATION`. The manual guide supplies the guided 3D central-brain ROI. Detection is performed on individual Z planes and detections are linked across Z. Missing manual annotations are reported as skipped samples.
4. Enable `RUN_QC` to generate QC from the saved detections. Raw IMS images are needed for image overlays.
5. Enable `RUN_SUMMARY` to export per-brain and group summaries and the configured comparisons. After a kernel restart, `existing_result_dirs` can select completed per-brain result directories for QC and reporting.

The maximum projection is an annotation aid; droplet counts come from the 3D workflow. The statistical unit in group summaries is a **brain**. Group definitions and the selected statistical method are stored in `workflow.json`.

Quantification creates new folders under `results/batches/`; reports go under `results/reports/`, unless another results root is configured. Per-brain outputs include `per_brain.csv`, `per_object.csv`, `per_plane.csv`, `per_detection.csv`, `config_used.yaml` and `run_info.json`. QC and marker stacks depend on the selected export options. Preserve the different definitions of LD-positive core volume and equivalent-sphere volume when interpreting the exported measurements.

## RNAscope quantification

This module includes **LD4 single-channel** and **LP dual-channel** analyses. Open [RNAscope_quantification_usage.ipynb](RNAscope%20analysis/RNAscope_quantification_usage.ipynb), choose `EXPERIMENT = 'LD4'` or `'LP_dual'`, and set `DATA_ROOT`.

### Configure datasets and masks

Populate `datasets` in the selected [LD4 configuration](RNAscope%20analysis/config/LD4.json) or [LP configuration](RNAscope%20analysis/config/LP_dual.json). The distributed lists are empty. Each entry has a unique ID and image/mask directories, for example:

```json
{
  "id": "batch01",
  "ims_dir": "LD4/batch01/ims",
  "mask_dir": "LD4/batch01/masks"
}
```

Relative dataset directories are resolved under `DATA_ROOT`. Each mask is paired with an IMS file having the same complete basename before the mask suffix:

| Analysis | Image channels, zero-based | Mask suffix |
| --- | --- | --- |
| LD4 | Channel 0 | `_CH1_mask.nrrd` |
| LP dual | Channels 0 and 1 | `_CH3_mask.nrrd`, or `_CH2_mask.nrrd` if no CH3 masks exist in that dataset |

For example, `sample01_LD4_ZT12.ims` pairs with `sample01_LD4_ZT12_CH1_mask.nrrd`. LD4 filenames must contain `ZT<number>`. LP filenames must contain the pattern `_LP_ZT<start>_<end>_<condition>_`, where the condition is a duration such as `60min` or a dark-control time such as `ZT16`.

Masks are externally annotated 3D integer label volumes in the same Z,Y,X geometry as the IMS images. Labels **1–4** identify sLNv cells, and label **5** identifies the 5th sLNv. The original counting rule assigns zero to absent labels; each sample contributes the mean count of labels 1–4 and the count of label 5.

### Quantification and outputs

Inspect `plan_experiment(config)` before enabling `RUN_QUANTIFICATION`. The supplied settings use 300-nm spot radii, Z,Y,X voxel spacing of 500 × 96 × 96 nm, and dense-region decomposition. LD4 uses threshold 100; LP uses channel-0 threshold 100 and channel-1 threshold 30. Match these parameters to the acquisition when analyzing other images.

Optional preprocessing exports a normalized structural channel for manual mask annotation. Spot detection reads the raw IMS channels. Optional marker export uses the coordinates saved during detection.

Each execution writes a new report directory under `results/<experiment>/`. Outputs include the XLSX report, plots when enabled, `configuration.json`, `input_manifest.csv`, `analysis.log` and `failures.csv`. Raw counts, analysis counts, excluded samples and statistics occupy separate report sheets. Configure exclusions with a regular-expression `pattern` and a `reason`; the supplied exclusion lists are empty.

The cache retains matching detections and spot-coordinate archives. Keep it together with the reports because `Spot_Files` references those archives by relative paths. Checkpoint reuse depends on input signatures, parameters, software and code. `random_seed = null` retains the original unspecified random state; an integer enables stable per-sample seeds for new analyses.

## Live-imaging cell tracking

The reference workflow is the locally used **Win10 adaptation of 3DeeCellTracker 0.5.2-alpha**. It uses the original method by Wen and colleagues, with replacement source files and models for the recordings analyzed here. Full installation details and the source comparison are in [Cell_tracking/README.md](Cell_tracking/README.md).

### Reconstruct the upstream installation

Use the working Win10 Conda environment; the source notebook records Python **3.8.19**. A complete historical TensorFlow/Keras/CUDA environment lock is not included. Consult the [pinned upstream installation instructions](https://github.com/WenChentao/3DeeCellTracker/blob/6819c761d586a1f66cd5615bcfd7fa9bcf5992f8/README.md) when preparing dependencies.

Activate that environment, then run the following from this repository's root in PowerShell. The checkout destination must be a new directory:

```powershell
$bundleDir = (Resolve-Path 'Cell_tracking').Path
$trackerSource = Join-Path $env:USERPROFILE 'source\3DeeCellTracker_win10'
New-Item -ItemType Directory -Path (Split-Path -Parent $trackerSource) -Force | Out-Null
git clone https://github.com/WenChentao/3DeeCellTracker.git "$trackerSource"
git -C "$trackerSource" checkout 6819c761d586a1f66cd5615bcfd7fa9bcf5992f8
Copy-Item -Path (Join-Path $bundleDir 'replacement_files\CellTracker\*.py') -Destination (Join-Path $trackerSource 'CellTracker') -Force
python -m pip install --no-deps -e "$trackerSource"
python -c "import CellTracker; print(CellTracker.__file__)"
```

`--no-deps` retains the dependencies in the existing working environment. The imported package should point to this patched checkout. Register/select that environment's kernel, and use **classic Jupyter Notebook** for the `%matplotlib notebook` tracking display. Restart the kernel after replacing package files. The notebook checks the replacement files and initial models before analysis.

### Select models and track a recording

| `MODEL_PROFILE` | Complete image XY size | Initial U-Net |
| --- | --- | --- |
| `128x128` | 128 × 128 pixels | `models/128x128/unet3_pretrained.h5` |
| `256x256` | 256 × 256 pixels | `models/256x256/unet3_pretrained.h5` |

Both profiles use `models/ffn_k3.hdf5`. The profiles describe complete image sizes; the serialized U-Net input patches are 80 × 80 × 5 and 80 × 80 × 16 voxels, respectively. Matching requires at least **four detected cells** in each point set, corresponding to the cell itself and three neighbors. `MIN_SIZE` is a region-size threshold in voxels.

1. Set `BUNDLE_DIR`, a separate `PROJECT_DIR` for the recording, `MODEL_PROFILE`, image naming, time/Z counts and voxel geometry in the notebook.
2. Run the setup sections and place one 2D TIFF per time/Z plane in `PROJECT_DIR/data/`. Indices start at 1; the default filename pattern is `Den_t%03i_z%03i.tif`.
3. Segment the first volume. Correct its labels in an external annotation tool and export one integer label TIFF per Z slice to `manual_vol1/`, ordered by filename. Label 0 is background, and each cell has a distinct positive integer ID consistent across Z slices.
4. Load the corrected first volume, fine-tune the selected initial U-Net, and select an actually saved checkpoint with `SELECTED_UNET_STEP` in the selection cell. The historical example used step 9; available checkpoints depend on the new training run.
5. Inspect the test matching, then track the remaining volumes and review the outputs.

Tracked labels are written to `track_results_SingleMode/`. Other outputs include `track_information/`, `anim/`, U-Net checkpoints and `run_config.json`. Preserve the corrected first-volume labels and the selected recording-specific weights with the recording. The bundled U-Nets are initial models for per-recording fine-tuning.

## Live-imaging MS2 signal processing

[Signal_processing.m](Signal_processing/Signal_processing.m) is a staged MATLAB entry function. The detailed review controls and output definitions are in [Signal_processing/README.md](Signal_processing/README.md).

### Prepare and configure a recording

Use one recording directory containing the prepared images and tracked masks:

```text
recording_01/
  data/                         # Input TIFF slices for Cell_tracking
  raw_data/                     # Input TIFF slices for Signal_processing
  manual_vol1/                  # Corrected first-volume labels
  models/                       # Recording-specific tracking weights
  track_results_SingleMode/     # Tracked label TIFF slices
```

Signal extraction requires single-plane **16-bit grayscale TIFF** images and integer label TIFF masks with matching dimensions. Names must end in `_t<number>_z<number>.tif` or `.tiff`, with complete time/Z grids starting at 1. Mask IDs must be consecutive `1:N`, with the same IDs present in every time volume; 0 denotes background.

Configure the acquisition once, using the actual values for the recording:

```matlab
repoDir = 'C:/path/to/this/repository';
addpath(fullfile(repoDir, 'Signal_processing'));
projectDir = 'C:/recordings/recording_01';

options = struct();
options.tvalue = 153;              % Example number of time points
options.zvalue = 43;               % Example Z slices per time point
options.xvalue = 256;              % Image rows
options.yvalue = 256;              % Image columns
options.StartTime = '4:07:00';      % Starting ZT/CT time
options.Cycle = 'LD3DD1';           % Recording label
options.FrameInterval = 44.42;     % Seconds between time points
options.signalsize = 5;            % Signal-window side length, pixels
options.bgsize = 20;               % Background-window side length, pixels
options.ci = 0.9;                  % Initial ON/OFF classification quantile
options.normalizationCells = [];  % Empty selects all cells
options.smooth_span = 0.05;        % Display LOESS span
options.plotOptimized = true;     % Include the five-method plots

Signal_processing('configure', projectDir, options);
Signal_processing('status', projectDir);
```

Create the recording directory before configuration. Configuration is saved in `SigInfo.mat`; subsequent stages use these saved settings.

If starting from ND2, run `Signal_processing('denoise', projectDir)` first. The default ND2 input is the sibling file `recording_01.nd2`; set `options.nd2File` during configuration to use another path. The reader follows the original two-channel acquisition layout. Prepared TIFF inputs allow this stage to be skipped.

For the tracking handoff, copy the prepared images from `raw_data/` to the tracker's `data/` directory. Set the tracking notebook to the same recording directory, matching time/Z counts, and `IMAGE_NAME = 'DenRaw_t%04i_z%04i.tif'` for MATLAB-denoised files. Retain `raw_data/` for MATLAB signal extraction. The tracker writes its labels directly into the recording's `track_results_SingleMode/` folder.

### Run stages individually

Execute each command separately and complete the manual steps before continuing:

| Command | Action |
| --- | --- |
| `Signal_processing('extract', projectDir)` | Extract per-cell image stacks and initial signal measurements |
| `Signal_processing('identify', projectDir)` | Assign cell types in the GUI and click Save |
| `Signal_processing('proofread', projectDir)` | Review positions and ON/OFF flags; save each cell before closing its window |
| `Signal_processing('correct', projectDir)` | Optional OFF-frame measurement correction, before normalization |
| `Signal_processing('normalize', projectDir)` | Apply baseline correction and normalization to the configured cells |
| `Signal_processing('optimize', projectDir)` | Compute the five projection/slice-selection methods |
| `Signal_processing('plot', projectDir)` | Export figures and numerical traces |
| `Signal_processing('ppt', projectDir)` | Optional PPT export; requires Report Generator and the complete five-method plot set |

Baseline correction uses robust LOESS with span 0.2; `smooth_span` controls display smoothing separately. If omitting the five-method branch, set `plotOptimized = false` during configuration and skip `optimize` and `ppt`.

Outputs remain in the recording directory, including `Cell_ID.mat`, `tracked_signal/`, `CellMaxIP/`, five-method outputs, and `SignalFigs/` containing figures and numerical MAT exports. Use `status` to inspect saved configuration and completed stages. To restart an analysis, use a fresh recording directory containing copies of the prepared inputs.

## Data and reproducibility

The figure spreadsheets are the inputs to the distributed MATLAB analyses. Raw IMS/ND2 images, recording-specific cell labels, and final fine-tuned tracking weights are needed to repeat the image-processing stages. The image-analysis configurations reference user-supplied data locations; they do not include a runnable raw-image example dataset.

Keep saved configurations, manual annotations, selected models and numerical outputs with each analysis. Bootstrap entry scripts provide their random seeds. Tracking fine-tuning and RNAscope decomposition retain the original unspecified seeds by default, with optional seed settings for new runs; numerical identity across different GPU/software environments is not guaranteed. Image workflows require execution and review with the matching images in their target environments.

## Citation and third-party components

Please cite the accompanying study when using these analyses. For cell tracking, also cite the original 3DeeCellTracker method and its U-Net foundation:

- Wen C, Miura T, Voleti V, et al. (2021). **3DeeCellTracker, a deep learning-based pipeline for segmenting and tracking cells in 3D time lapse images.** *eLife* 10:e59187. https://doi.org/10.7554/eLife.59187
- Çiçek Ö, Abdulkadir A, Lienkamp SS, Brox T, Ronneberger O. (2016). **3D U-Net: Learning Dense Volumetric Segmentation from Sparse Annotation.** MICCAI. https://doi.org/10.1007/978-3-319-46723-8_49

The tracking code is based on [WenChentao/3DeeCellTracker](https://github.com/WenChentao/3DeeCellTracker), with historical adaptation by [vislabwwy/3DCT_for_PKU_LDG](https://github.com/vislabwwy/3DCT_for_PKU_LDG) and the supplied local replacement files. The exact upstream reconstruction uses commit `6819c761d586a1f66cd5615bcfd7fa9bcf5992f8` as specified above. RNAscope detection uses [Big-FISH](https://github.com/fish-quant/big-fish). The bundled `nd2finfo.m` retains its original Chao-Yuan Yeh copyright notice.

Third-party components retain their own license and copyright notices. The upstream 3DeeCellTracker MIT license is preserved in [Cell_tracking/LICENSE_3DeeCellTracker](Cell_tracking/LICENSE_3DeeCellTracker). The supplied model files did not include a separate model-license statement; model reuse terms should be confirmed with their providers before redistribution. A repository-level license does not replace those notices or establish separate model reuse rights.
