function SigInfo = Signal_processing(stage, projectDir, options)
%SIGNAL_PROCESSING Run one stage of the live-imaging MS2 signal workflow.
% No arguments: show help only. Edit defaultConfiguration below, or pass
% configuration overrides in the third argument of the configure stage.
%
% Example (replace the path and acquisition parameters):
%   p = 'D:/recordings/example';
%   Signal_processing('configure',p,struct('tvalue',153,'zvalue',43));
%   Signal_processing('denoise',p); % optional if raw_data already exists
%   % Run Cell_tracking; provide label TIFFs in track_results_SingleMode.
%   Signal_processing('extract',p);
%   Signal_processing('identify',p); % select cell types, then click Save
%   Signal_processing('proofread',p);% Save each cell before closing its GUI
%   Signal_processing('correct',p);  % optional: recalculate OFF frames
%   Signal_processing('normalize',p);
%   Signal_processing('optimize',p); % five original projection methods
%   Signal_processing('plot',p);     % optional; includes numerical exports
%   Signal_processing('ppt',p);      % optional; MATLAB Report Generator
%
% Review selected cells/frames:
%   Signal_processing('proofread',p,struct('cells',[2 3],'frames',101:153));
% 'status' shows saved configuration and completed stages. Each stage is a
% separate call because identification and proofreading are manual. Existing
% extraction/normalization outputs are guarded; there is no run-all stage.
% Use a clean analysis directory when starting from the raw images.
%
% Requirements: MATLAB, Image Processing Toolbox, Statistics and Machine
% Learning Toolbox, Curve Fitting Toolbox (smooth/rloess/loess). PPT additionally
% needs MATLAB Report Generator. See README.md for setup and usage.

if nargin<1 || isempty(stage), stage='help'; end
if nargin<2, projectDir=''; end
if nargin<3, options=struct(); end
stage=validatestring(stage,{'help','configure','status','denoise','extract', ...
    'identify','proofread','correct','normalize','optimize','plot','ppt'});
SigInfo=struct();
if strcmp(stage,'help'), help(mfilename); return; end
validateattributes(options,{'struct'},{'scalar'});
if isempty(projectDir)
    projectDir=uigetdir(pwd,'Select the recording directory');
    if isequal(projectDir,0), return; end
end
projectDir=char(projectDir);
assert(isfolder(projectDir),'MS2:MissingProject','Missing recording directory: %s',projectDir);
[~,attr]=fileattrib(projectDir); projectDir=attr.Name;
addpath(fullfile(fileparts(mfilename('fullpath')),'script'));
previousDir=pwd;
restoreDir=onCleanup(@() cd(previousDir)); %#ok<NASGU>
cd(projectDir); % retained functions use relative output paths
configFile=fullfile(projectDir,'SigInfo.mat');
if strcmp(stage,'configure')
    assert(~isfile(configFile),'MS2:ExistingConfiguration', ...
        'SigInfo.mat exists. Inspect status; do not overwrite an existing analysis.');
    SigInfo=defaultConfiguration();
    keys=fieldnames(options);
    for k=1:numel(keys)
        assert(isfield(SigInfo,keys{k}),'MS2:UnknownOption','Unknown option: %s',keys{k});
        SigInfo.(keys{k})=options.(keys{k});
    end
    validateConfiguration(SigInfo);
    SigInfo.p=projectDir;
    SigInfo.PublicWorkflowVersion=1;
    SigInfo.CompletedStages={'configure'};
    save(configFile,'SigInfo');
    fprintf('Configuration saved in %s. Inspect it with status.\n',configFile);
    return
end
assert(isfile(configFile),'MS2:ConfigureFirst','Run configure first.');
saved=load(configFile,'SigInfo'); SigInfo=saved.SigInfo;
assert(isfield(SigInfo,'PublicWorkflowVersion'),'MS2:LegacyConfiguration', ...
    'Legacy SigInfo.mat: configure this workflow in a separate clean analysis directory.');
SigInfo.p=projectDir; % allow moving a complete recording directory
if ~strcmp(stage,'proofread')
    assert(isempty(fieldnames(options)),'MS2:FrozenConfiguration', ...
        'Pass acquisition overrides at configure; only proofread accepts per-run cells/frames.');
end
if strcmp(stage,'status'), disp(SigInfo); return; end
switch stage
    case 'denoise'
        requireEmpty(fullfile(projectDir,SigInfo.rawfolder));
        DeNoise(SigInfo);
    case 'extract'
        for field={'targetfolder','signalfolder','MIPfolder'}
            requireEmpty(fullfile(projectDir,SigInfo.(field{1})));
        end
        assert(~isfile('Cell_ID.mat') && ~isfile('cell_labels.tif'), ...
            'MS2:ExistingExtraction','Cell identification/extraction outputs already exist.');
        SigInfo.totalcellnum=TrackedDataSave(SigInfo);
        for i=1:SigInfo.totalcellnum
            SignalExtraction_v9(SigInfo,i,SigInfo.signalsize,SigInfo.bgsize,SigInfo.ci);
        end
    case 'identify'
        requireStage(SigInfo,'extract');
        if ~isfile('Cell_ID.mat'), fig=Cellid(SigInfo.totalcellnum); waitfor(fig); end
        assert(isfile('Cell_ID.mat'),'MS2:IdentificationNotSaved', ...
            'Identification closed without saving; run identify again.');
        renamefolder('Cell_ID.mat',SigInfo.signalfolder);
    case 'proofread'
        requireStage(SigInfo,'identify');
        cells=SigInfo.proofreadCells; frames=SigInfo.proofreadFrames;
        assert(all(ismember(fieldnames(options),{'cells','frames'})), ...
            'MS2:UnknownOption','Proofread accepts cells and frames only.');
        if isfield(options,'cells'), cells=options.cells; end
        if isfield(options,'frames'), frames=options.frames; end
        cells=resolveCells(cells,SigInfo.totalcellnum);
        if isempty(frames), frames=1:SigInfo.tvalue; end
        validateattributes(frames,{'numeric'},{'row','integer','positive','increasing','<=',SigInfo.tvalue});
        assert(all(diff(frames)==1),'MS2:FrameRange','Use a continuous frame range.');
        [~,~,files]=listSignalFolders(fullfile(projectDir,SigInfo.signalfolder),cells);
        requireUnnormalized(files);
        fprintf('For each cell, click Save before closing its proofreading window.\n');
        for k=1:numel(files)
            data=load(files{k},'CellSignal');
            data.CellSignal.path=projectDir; data.CellSignal.SigInfo.p=projectDir;
            fig=GenFigure_MIPsig_ci(data.CellSignal,frames); waitfor(fig);
        end
        fprintf('Review windows closed. Proceed after saving all intended corrections.\n');
        return % opening a GUI does not certify completion of manual review
    case 'correct'
        requireStage(SigInfo,'identify'); requireNotCompleted(SigInfo,'correct');
        [~,~,files]=listSignalFolders(fullfile(projectDir,SigInfo.signalfolder));
        requireUnnormalized(files);
        for i=1:SigInfo.totalcellnum, furtherCorrection_v3(SigInfo,i); end
    case 'normalize'
        requireStage(SigInfo,'identify'); requireNotCompleted(SigInfo,'normalize');
        cells=resolveCells(SigInfo.normalizationCells,SigInfo.totalcellnum);
        [~,~,files]=listSignalFolders(fullfile(projectDir,SigInfo.signalfolder),cells);
        requireUnnormalized(files);
        for i=cells, NormalizeResults_spline250822_all_cell_min_20x20(SigInfo,i); end
    case 'optimize'
        requireStage(SigInfo,'normalize'); requireNotCompleted(SigInfo,'optimize');
        requireEmpty(fullfile(projectDir,'ProcessedBySF_traces'));
        assert(~isfile('Hash_ymaxV.mat'),'MS2:ExistingOptimization','Hash_ymaxV.mat already exists.');
        OptimizeResults(SigInfo);
    case 'plot'
        requireStage(SigInfo,'normalize');
        if SigInfo.plotOptimized, requireStage(SigInfo,'optimize'); end
        % Legacy plots need the post-normalization eight-column MIPsig layout.
        plotInfo=SigInfo;
        plotInfo.plotCells=resolveCells(SigInfo.normalizationCells,SigInfo.totalcellnum);
        PlotSignals(plotInfo); PlotSignals_raw(plotInfo); PlotNormalizedSignals(plotInfo);
        PlotSingleCellFig(plotInfo); PlotSingleCellFig_raw(plotInfo);
        PlotNormalizedSingleCellFig(plotInfo);
        PlotNormalizedSingleCellFig_all20x20(plotInfo);
        PlotNormalizedSingleCellFig_partial20x20(plotInfo);
        PlotSignals_loess(plotInfo,SigInfo.smooth_span);
        PlotNormalizedSignals_loess(plotInfo,SigInfo.smooth_span);
        PlotSingleCellFig_loess(plotInfo,SigInfo.smooth_span);
        PlotNormalizedSingleCellFig_loess(plotInfo,SigInfo.smooth_span);
        if SigInfo.plotOptimized
            PlotOptimizedSignals(plotInfo); PlotOptimizedSingleCellFig(plotInfo);
            PlotOptimizedSignals_loess(plotInfo,SigInfo.smooth_span);
            PlotOptimizedSingleCellFig_loess(plotInfo,SigInfo.smooth_span);
        end
    case 'ppt'
        requireStage(SigInfo,'plot');
        assert(SigInfo.plotOptimized,'MS2:IncompletePPT','The PPT template requires the optimized plots.');
        [~,name]=fileparts(projectDir);
        assert(~isfile([name '.pptx']),'MS2:ExistingPPT','The recording PPT already exists.');
        saveppt_loess_20x20(projectDir);
end
SigInfo.CompletedStages=unique([SigInfo.CompletedStages,{stage}],'stable');
save(configFile,'SigInfo');
fprintf('Stage %s finished.\n',stage);
end

function cfg=defaultConfiguration()
% Acquisition values are examples from the original entry, not inferred.
cfg.tvalue=153; cfg.zvalue=43;
cfg.xvalue=256; cfg.yvalue=256; % historical names: image rows, image columns
cfg.StartTime='4:07:00'; cfg.Cycle='LD3DD1';
cfg.FrameInterval=44.42; % SECONDS
cfg.rawfolder='raw_data'; cfg.trackfolder='track_results_SingleMode';
cfg.targetfolder='tracked_raw_data'; cfg.signalfolder='tracked_signal';
cfg.spotfolder='detected_spot'; cfg.MIPfolder='CellMaxIP';
cfg.nd2File=''; % empty: sibling <recording directory name>.nd2
cfg.signalsize=5; cfg.bgsize=20; cfg.ci=0.9;
cfg.normalizationCells=[]; % empty: all; old entry selected [2 3 5 7]
cfg.proofreadCells=[];     % empty: all
cfg.proofreadFrames=[];    % empty: all; old launchers used 1:100
cfg.smooth_span=0.05;      % display LOESS; baseline rloess stays 0.2
cfg.plotOptimized=true;   % false: omit five-method plots and PPT
end

function validateConfiguration(cfg)
for field={'tvalue','zvalue','xvalue','yvalue','signalsize','bgsize'}
    validateattributes(cfg.(field{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
validateattributes(cfg.FrameInterval,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.ci,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(cfg.smooth_span,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.plotOptimized,{'logical'},{'scalar'});
datevec(cfg.StartTime);
folders={'rawfolder','trackfolder','targetfolder','signalfolder','spotfolder','MIPfolder'};
names=cell(size(folders));
for k=1:numel(folders)
    names{k}=char(cfg.(folders{k}));
    assert(~isempty(regexp(names{k},'^[A-Za-z0-9_ -]+$','once')), ...
        'MS2:FolderName','Use a single child-directory name for %s.',folders{k});
end
assert(numel(unique(lower(string(names))))==numel(names),'MS2:FolderName','Directory names must be distinct.');
for field={'normalizationCells','proofreadCells','proofreadFrames'}
    value=cfg.(field{1});
    if ~isempty(value)
        validateattributes(value,{'numeric'},{'row','integer','positive','finite'});
        assert(numel(unique(value))==numel(value),'MS2:DuplicateSelection','Duplicate selection.');
    end
end
end

function cells=resolveCells(cells,total)
if isempty(cells), cells=1:total; end
validateattributes(cells,{'numeric'},{'row','integer','positive','<=',total});
assert(numel(unique(cells))==numel(cells),'MS2:DuplicateSelection','Duplicate cell selection.');
end

function requireEmpty(folder)
if isfolder(folder)
    entries=dir(folder);
    assert(all(ismember({entries.name},{'.','..'})),'MS2:ExistingOutput', ...
        'Output is not empty: %s. Use a clean analysis directory or move prior outputs aside.',folder);
end
end

function requireStage(info,stage)
assert(ismember(stage,info.CompletedStages),'MS2:StageOrder','Run stage %s first.',stage);
end

function requireNotCompleted(info,stage)
assert(~ismember(stage,info.CompletedStages),'MS2:RepeatedStage', ...
    'Stage %s already finished; numerical outputs are not overwritten.',stage);
end

function requireUnnormalized(files)
for k=1:numel(files)
    data=load(files{k},'CellSignal');
    assert(~isfield(data.CellSignal,'NormalizedMIPsig') && size(data.CellSignal.MIPsig,2)==7, ...
        'MS2:NormalizedInput','Expected a seven-column pre-normalization result: %s',files{k});
end
end
