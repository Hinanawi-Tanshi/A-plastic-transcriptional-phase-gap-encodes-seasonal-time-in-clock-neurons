# Cell tracking

This folder packages the laboratory's **Win10, U-Net-based single-mode workflow** for segmenting and tracking cells in live 3D time-lapse recordings. It follows the locally used `Examples/use_unet/single_mode_worm1-clear.ipynb` and supplies the files needed to reconstruct its customized 3DeeCellTracker installation.

The workflow is: initial segmentation of the first 3D volume, manual correction of that volume, U-Net fine-tuning using the corrected labels, then FFN + PR-GLS tracking of subsequent volumes. Each recording needs its own images and first-volume labels.

## Contents

```text
Cell_tracking/
  Cell_tracking.ipynb
  README.md
  LICENSE_3DeeCellTracker
  replacement_files/
    CellTracker/
      tracker.py
      unet3d.py
      ffn.py
      synthesize.py
  models/
    128x128/unet3_pretrained.h5
    256x256/unet3_pretrained.h5
    ffn_k3.hdf5
```

The replacement files are unchanged copies of the locally used Win10 sources. The notebook has been reorganized around one recording configuration, with saved outputs cleared. It checks the installed replacement files and bundled models using SHA-256 checksums before starting the analysis.

## 1. Prepare the environment

Use the **working Win10 tracking Conda environment** used for the original analysis. The source notebook records Python **3.8.19**. A complete export of that machine's TensorFlow/Keras/CUDA environment was not supplied, so this folder does not provide a verified fresh-machine environment lock. The Win11 environment is not used as the reference environment here.

The upstream prerequisites are TensorFlow 2.x, its compatible CUDA/cuDNN installation for GPU use, and the scientific Python/Jupyter dependencies of 3DeeCellTracker. For a new installation, consult the [README at the pinned upstream commit](https://github.com/WenChentao/3DeeCellTracker/blob/6819c761d586a1f66cd5615bcfd7fa9bcf5992f8/README.md). Keep the already working environment when reconstructing an existing analysis. The upstream `requirements.txt` is not a verified lock for the historical Win10 machine.

Use **classic Jupyter Notebook** for the original `%matplotlib notebook` interactive tracking display. Activate the intended environment in an Anaconda PowerShell Prompt before following the installation commands below. Do not run these commands from a different Python environment.

## 2. Get the exact upstream source and replace four files

The reference is [WenChentao/3DeeCellTracker](https://github.com/WenChentao/3DeeCellTracker), **0.5.2-alpha**, at commit:

```text
6819c761d586a1f66cd5615bcfd7fa9bcf5992f8
```

Use this commit rather than the current default branch. The replacement files were compared against this exact snapshot. Keep the source checkout separate from this publication folder.

From the activated Win10 environment, run the following in PowerShell, changing the first path to this folder's location:

```powershell
Set-Location 'C:\path\to\Cell_tracking'
$bundleDir = (Get-Location).Path
$trackerSource = Join-Path $env:USERPROFILE 'source\3DeeCellTracker_win10'
$upstreamCommit = '6819c761d586a1f66cd5615bcfd7fa9bcf5992f8'

New-Item -ItemType Directory -Path (Split-Path -Parent $trackerSource) -Force | Out-Null
git clone https://github.com/WenChentao/3DeeCellTracker.git "$trackerSource"
git -C "$trackerSource" checkout $upstreamCommit

Copy-Item -Path (Join-Path $bundleDir 'replacement_files\CellTracker\*.py') `
    -Destination (Join-Path $trackerSource 'CellTracker') -Force

python -m pip install --no-deps -e "$trackerSource"
python -c "import CellTracker; print(CellTracker.__file__)"
```

The clone destination must be a new directory. `--no-deps` retains the dependencies in the working environment; it assumes those dependencies are already installed. The printed package path should point into this patched source checkout.

Register the environment as a notebook kernel if needed:

```powershell
python -m ipykernel install --user --name cell_tracking_win10 --display-name 'Cell tracking (Win10)'
python -m notebook Cell_tracking.ipynb
```

Select **Cell tracking (Win10)**. Restart the kernel after replacing source files. The notebook stops if Python imports a package whose four replacement files differ from the supplied copies.

### Changes relative to the pinned upstream source

| Replacement file | Preserved Win10 change | Role in this workflow |
|---|---|---|
| `tracker.py` | FFN matching neighbor count changed from 20 to **3** | Active tracking path; use the supplied three-neighbor FFN |
| `unet3d.py` | `unet3_a()` factory input changed from `(160,160,16,1)` to `(32,160,16,1)` | Historical factory change; this notebook loads the supplied serialized U-Net models, whose input shapes are listed below |
| `ffn.py` | Training-side `k_ptrs` changed from 20 to 13 | Preserved historical source; this notebook loads the supplied FFN and does not train one |
| `synthesize.py` | Training-side `K_NEIGHBORS` changed from 20 to 13 | Preserved historical source; not the neighbor setting used by `Tracker` during tracking |

The training-side constants are not a ready-to-use recipe for retraining the three-neighbor FFN. They are retained to reconstruct the historical code. The tracking path uses `CellTracker.track.initial_matching_quick(..., 3)` and the supplied 10-feature FFN.

## 3. Select the models

| Notebook profile | Complete image XY size | Serialized U-Net input patch, excluding batch/channel | Bundled U-Net |
|---|---|---|---|
| `128x128` | 128 x 128 pixels | 80 x 80 x 5 voxels | `models/128x128/unet3_pretrained.h5` |
| `256x256` | 256 x 256 pixels | 80 x 80 x 16 voxels | `models/256x256/unet3_pretrained.h5` |

Both profiles use **`models/ffn_k3.hdf5`**. This is the historical `weightsweights.26.hdf5` file, renamed to identify its role. The copies in the original 128 and 256 model folders were byte-identical. Each of its two input branches has 10 features: `3 neighbors x 3 coordinates + 1 distance feature`.

The other file named `ffn_pretrained.h5` in the historical model folders has 16-feature inputs and is not used by these notebooks; it is not included here. Intermediate retraining checkpoints and example-specific `manual_vol1.nrrd` files are also not bundled as general initial models.

The notebook copies the selected initial U-Net and the common FFN into the recording's own `models` directory. Existing files are reused only when their checksums match. It fine-tunes the loaded U-Net directly; it does not call the legacy `unet3_a()` factory.

**Cell-count condition:** both point sets used in matching must contain at least **four detected cells** because the query includes the cell itself and its three neighbors. This is an input-count requirement, not a guarantee of tracking accuracy. `min_size=100` is a region-size threshold in voxels, not a minimum number of cells.

## 4. Configure and prepare a recording

Edit section 1 of `Cell_tracking.ipynb`:

- Set `BUNDLE_DIR` to the folder containing the notebook if it is not the notebook's working directory.
- Choose `MODEL_PROFILE` once, and use a distinct `PROJECT_DIR` for each recording.
- Set `VOLUME_NUM`, `Z_SLICES` and `IMAGE_NAME` for the supplied images.
- Set `Z_XY_RATIO`, `Z_SCALING`, `NOISE_LEVEL`, `MIN_SIZE` and the tracking parameters for the recording.

The default 256-mode settings follow the source Win10 notebook: 135 volumes, 46 Z slices, `z_xy_ratio=4.98`, `z_scaling=5`, noise level 20 and minimum region size 100. These are editable example parameters. Selecting the 128 profile changes XY size and U-Net selection; update the Z depth and all recording-specific settings as needed.

Run sections 1-3 to verify the patched package, select the models and create the project directories. Copy the images into the created `PROJECT_DIR/data` directory:

```text
experiment_01/
  data/
    Den_t001_z001.tif
    Den_t001_z002.tif
    ...
    Den_t135_z046.tif
  manual_vol1/
  models/
```

Each TIFF is one 2D slice. Time and Z numbering start at **1**. Match the filename pattern, complete XY dimensions and Z depth to the notebook settings. Section 4 checks the expected files and segments the first volume.

## 5. Correct the first volume and fine-tune U-Net

1. Inspect the initial segmentation. In an external annotation tool such as ITK-SNAP, correct the first 3D volume: remove false detections, add missed cells and refine boundaries.
2. Export the corrected label images into `PROJECT_DIR/manual_vol1`: one 2D label TIFF per Z slice, sorted in Z order by filename, such as `label_z001.tif`. Use **0 for background and a distinct positive integer for each cell**, preserving the same ID across Z slices. Export label values rather than RGB overlays. Keep only these label TIFFs in that directory.
3. Run section 5 to validate and load the manual labels. Their dimensions must match the first imaging volume, and the reference must contain at least four cells.
4. Run section 6 to fine-tune the selected initial U-Net on this recording's corrected first volume. The original routine uses 10 training steps and saves checkpoints only when validation loss improves.
5. Review the predictions and available checkpoints. Set `SELECTED_UNET_STEP` in the **selection cell in section 6** and run that cell. The historical Win10 example selected step **9**, but that step may not be saved during a new training run. Choose a checkpoint actually produced for the current recording. Step 0 selects the initial weights.

When changing `SELECTED_UNET_STEP`, rerun the selection cell directly; rerunning the configuration cell is unnecessary. To change the input recording or initial model profile, restart the kernel and use a fresh project directory.

The training cell stops if fine-tuning checkpoints already exist, so stale checkpoints from another training run are not selected accidentally. Use a fresh project folder for a new training run. To reuse a deliberately archived checkpoint for the same recording and manual labels, skip the training cell and select that checkpoint explicitly in the following cell.

The selection cell clears the U-Net prediction cache, as in the original notebook, so subsequent predictions use the selected weights. Training was not seeded in the historical notebook. `RANDOM_SEED=None` preserves that behavior; optionally set a seed for a new run. GPU training need not be numerically identical across environments even with a seed.

Preserve the first-volume labels and **the final selected recording-specific weights** with the analysis. The two bundled U-Nets are initial models, not a substitute for the weights selected after per-recording fine-tuning.

## 6. Test matching and track subsequent volumes

Run section 7 to interpolate the manual labels, initialize tracking and load the three-neighbor FFN. The source Win10 tracking settings are `beta_tk=300`, `lambda_tk=0.2`, `maxiter_tk=500`. Inspect the match to `MATCH_VOLUME` (default 7), its correction display and overlapping regions.

Then run section 8 to track volumes 2 through `VOLUME_NUM`. Use classic Notebook for the interactive Matplotlib display. Review the replay animation and tracking outputs.

If a frame produces `Expected n_neighbors <= n_samples` with `n_neighbors=4`, inspect its images and segmentation: fewer than four cells were detected in a point set used for matching. A reference mask containing four cells does not ensure every later frame will detect four. Check cell detection and tracking quality before accepting the result. Do not change the neighbor count independently of the FFN input dimensions.

## Outputs

All recording-specific outputs are written under `PROJECT_DIR` using the upstream directory layout:

| Directory/file | Content |
|---|---|
| `auto_vol1/` | Automatic first-volume segmentation |
| `manual_vol1/` | User-supplied first-volume labels |
| `models/unet_weights/` | Initial weights, saved fine-tuning checkpoints and selected model |
| `unet_cache/` | Cached U-Net predictions |
| `track_results_SingleMode/` | Tracked label slices |
| `track_information/` | Tracking coordinates/information |
| `anim/` | Tracking process images |
| `run_config.json` | Effective parameters, source/model/label hashes, selected checkpoint and run stage |

The distributed notebook has no saved execution output or image data. Raw imaging data, manual annotations and per-recording trained weights must be retained with the corresponding experiment.

## Provenance, citation and third-party license

- Original method and implementation: [WenChentao/3DeeCellTracker](https://github.com/WenChentao/3DeeCellTracker), pinned as specified above.
- Historical Luo-lab adaptation: [vislabwwy/3DCT_for_PKU_LDG](https://github.com/vislabwwy/3DCT_for_PKU_LDG). That public repository is based on version 0.4.5 and is credited as development history. The supplied replacement files reconstruct the later locally used Win10 version against the pinned 0.5.2-alpha source.
- The replacement sources retain the original author's notices. The upstream MIT license is preserved in `LICENSE_3DeeCellTracker`. The separately supplied model assets did not include a model-specific license statement; their redistribution and reuse terms should be recorded with their providers before external publication.

Please cite:

Wen C, Miura T, Voleti V, et al. (2021). **3DeeCellTracker, a deep learning-based pipeline for segmenting and tracking cells in 3D time lapse images.** *eLife* 10:e59187. https://doi.org/10.7554/eLife.59187

This workflow uses the U-Net route. The upstream project also recommends citing the 3D U-Net method: Cicek O, Abdulkadir A, Lienkamp SS, Brox T, Ronneberger O. (2016). **3D U-Net: Learning Dense Volumetric Segmentation from Sparse Annotation.** MICCAI. https://doi.org/10.1007/978-3-319-46723-8_49

## Validation of this public bundle

The preparation checks verify that the four replacement files and three model files match their historical source files, that applying the replacements reconstructs the local Win10 package relative to the pinned upstream snapshot, and that the notebook is valid and its code cells compile. Model input shapes are checked from HDF5 metadata.

These preparation checks do not constitute an end-to-end execution on the original Win10 imaging environment. Full GPU tracking and per-recording fine-tuning should be verified there with the matching image data and manual annotations.
