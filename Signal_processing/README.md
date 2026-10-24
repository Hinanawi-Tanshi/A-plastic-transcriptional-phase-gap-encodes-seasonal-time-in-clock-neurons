# MS2 signal processing

MATLAB workflow for extracting and classifying MS2 signals from tracked cells, manually reviewing detections, correcting the background baseline, and exporting traces and figures.

## Requirements

- MATLAB (tested with R2023a).
- Image Processing Toolbox.
- Statistics and Machine Learning Toolbox.
- Curve Fitting Toolbox.
- MATLAB Report Generator, only for PPT export.

`Signal_processing.m` is the entry function. It adds the accompanying `script` directory to the MATLAB path automatically. Run `Signal_processing` without arguments to display help.

## Prepare the recording directory

Keep each recording in its own directory, separate from the code:

```text
recordings/
├── example.nd2                         # optional input for denoise
└── example/                            # projectDir
    ├── raw_data/
    │   ├── DenRaw_t0001_z0001.tif
    │   └── ...
    └── track_results_SingleMode/
        ├── track_results_t0001_z0001.tif
        └── ...
```

The `extract` stage requires prepared images in `raw_data` and tracked masks from the companion `Cell_tracking` workflow in `track_results_SingleMode`.

- Images must be single-plane, 16-bit grayscale TIFFs. Masks must be single-plane integer label TIFFs with the same dimensions.
- Names must end in `_t<number>_z<number>.tif` or `.tiff`. Both `Den_t001_z001.tif` and `DenRaw_t0001_z0001.tif` are accepted.
- Time and Z indices start at 1. Both directories must contain the complete configured time-by-Z grid, with one file per plane.
- Mask label 0 is background; cell IDs are consecutive integers `1:N`. The same cell IDs must occur in every time volume.

If prepared images already exist, skip `denoise`; an ND2 file is then unnecessary. If using `denoise`, the default ND2 location is the sibling file `<recording directory name>.nd2`, as shown above. The bundled reader supports the original two-channel acquisition; check the channel assignment when using a different acquisition layout.

## Configure once per recording

Open this code directory in MATLAB, or add it to the path. Replace the example paths and acquisition values below with those for your recording:

```matlab
addpath('D:/path/to/Signal_processing');
projectDir = 'D:/recordings/example';

options = struct();
options.tvalue = 153;             % number of time points
options.zvalue = 43;              % Z slices per time point
options.xvalue = 256;             % image rows
options.yvalue = 256;             % image columns
options.StartTime = '4:07:00';     % starting ZT/CT time
options.Cycle = 'LD3DD1';          % recording label
options.FrameInterval = 44.42;    % seconds between time points
options.signalsize = 5;           % signal-window side length, pixels
options.bgsize = 20;              % background-window side length, pixels
options.ci = 0.9;                 % KDE quantile for initial ON/OFF detection
options.normalizationCells = []; % empty: all cells; e.g. [2 3 5 7]
options.smooth_span = 0.05;       % LOESS span for display
options.plotOptimized = true;    % include plots of the five-method branch

Signal_processing('configure', projectDir, options);
Signal_processing('status', projectDir);
```

Configuration is saved as `SigInfo.mat` in the recording directory. Additional defaults, including input/output folder names and `nd2File`, are in `defaultConfiguration` inside the entry function. To use an ND2 file elsewhere, include `options.nd2File = 'D:/recordings/input.nd2'` before calling `configure`.

Start with a recording directory containing the prepared inputs and no previous processing outputs or `SigInfo.mat`. Later stages use the saved configuration; editing the entry's defaults does not change an already configured recording.

## Run the workflow in stages

Execute each command separately. Finish and save the manual steps before continuing.

### 1. Prepare images and tracked masks

If starting from ND2, generate the denoised images:

```matlab
Signal_processing('denoise', projectDir);
```

Then run the companion `Cell_tracking` notebook. Set its `PROJECT_DIR` to this same recording directory, `IMAGE_NAME` to `DenRaw_t%04i_z%04i.tif`, and `VOLUME_NUM` and `Z_SLICES` to the MATLAB time and Z counts. Its tracked label images belong in `track_results_SingleMode`.

If images and masks are already prepared, proceed directly to extraction:

```matlab
Signal_processing('extract', projectDir);
```

### 2. Identify cells

```matlab
Signal_processing('identify', projectDir);
```

Select the cell type for each tracked ID and click **Save**. The workflow saves `Cell_ID.mat` and adds the cell type to each result folder's name. Repeating this stage with an unchanged mapping reuses the existing names.

### 3. Review detections manually

```matlab
Signal_processing('proofread', projectDir);
```

The workflow opens one cell's review window at a time. The frame buttons toggle ON/OFF status. For a position correction, select the annotated frame image and use **Correction**: a left click marks ON and a right click marks OFF. Changing **ci Thresh** and clicking **Update** recalculates classification over the displayed frame range. Click **Save** before closing the window.

By default, all cells and frames are shown. To review a subset:

```matlab
Signal_processing('proofread', projectDir, ...
    struct('cells', [2 3], 'frames', 101:153));
```

Frame ranges must be continuous and within the recording. Complete all intended review before normalization; the review interface operates on pre-normalization results.

### 4. Correct and normalize signals

Optionally recalculate OFF-frame measurements near the mean ON-frame signal position. This preserves the reviewed ON/OFF flags:

```matlab
Signal_processing('correct', projectDir);
```

Normalize the cells selected by `normalizationCells`:

```matlab
Signal_processing('normalize', projectDir);
```

Background baseline correction uses robust LOESS (`rloess`, span 0.2). The pre-normalization measurements are retained as `CellSignal.MIPsigBeforeNormalization`.

### 5. Run the five-method calculation branch

```matlab
Signal_processing('optimize', projectDir);
```

This computes the original maximum/mean projection and slice-selection methods for all extracted cells. To omit this branch, set `plotOptimized = false` during configuration and skip this command.

### 6. Export figures and an optional PPT

```matlab
Signal_processing('plot', projectDir);
```

Plots use the normalized cell selection. The display LOESS span is controlled by `smooth_span`, separately from background baseline correction. Numerical exports are also saved in `SignalFigs`.

To collect the figures into a presentation:

```matlab
Signal_processing('ppt', projectDir);
```

PPT export requires MATLAB Report Generator and the complete plot set, including the five-method plots (`plotOptimized = true`).

## Outputs

All outputs are written inside `projectDir`.

| Output | Contents |
| --- | --- |
| `SigInfo.mat` | Saved configuration and completed processing stages |
| `cell_labels.tif`, `Cell_ID.mat` | Tracked cell labels and selected cell types |
| `tracked_raw_data/` | Masked images for individual cells |
| `CellMaxIP/` | Per-cell maximum-intensity projection TIFF stacks |
| `tracked_signal/` | Per-cell `CellSignal` MAT files and review launchers |
| `ProcessedBySF_traces/`, `Hash_ymaxV.mat` | Five-method calculation outputs |
| `SignalFigs/` | FIG/PNG figures, `timetrace.mat`, and signal MAT exports |
| `<recording name>.pptx` | Optional presentation |

In `CellSignal.NormalizedMIPsig`, columns 1–5 contain the ON/OFF flag, baseline-corrected intensity, signal-normalized trace, background-referenced signal trace, and normalized background trace. The retained `NormalizedSignalFigures` overview displays the background column; its LOESS counterpart displays the signal-normalized column. `Allsignal.mat` contains full `CellSignal` structures, rather than only the plotted values.

## Reopening or repeating an analysis

Use `Signal_processing('status', projectDir)` to inspect the saved settings and stage history. Extraction, normalization, optimization, and PPT export guard against existing outputs. To restart processing, use a fresh recording directory with copies of the prepared images and masks.

Generated per-cell review launchers can also be opened independently before normalization. Add this bundle's `script` directory to the MATLAB path first. They accept an optional frame range and locate their MAT file relative to their own location.
