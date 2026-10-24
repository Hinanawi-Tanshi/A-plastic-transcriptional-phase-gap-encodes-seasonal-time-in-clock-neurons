function SaveMFile_MIPsig_ci(filename)
% Generate a movable per-cell launcher; no argument means all frames.
% Ensure the public script directory is on the MATLAB path before using it.
[folder,name]=fileparts(filename);
assert(isvarname(name),'MS2:LauncherName','Invalid function name: %s',name);
fid=fopen(fullfile(folder,[name '.m']),'w');
assert(fid>=0,'MS2:LauncherWrite','Cannot create proofreading launcher.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'function fig = %s(trange)\n',name);
fprintf(fid,'folder = fileparts(mfilename(''fullpath''));\n');
fprintf(fid,'files = dir(fullfile(folder,''cell_*_tracked_signal.mat''));\n');
fprintf(fid,'assert(numel(files)==1,''Expected one tracked-signal MAT file.'');\n');
fprintf(fid,'data = load(fullfile(files.folder,files.name),''CellSignal'');\n');
fprintf(fid,'data.CellSignal.path = fileparts(fileparts(folder));\n');
fprintf(fid,'data.CellSignal.SigInfo.p = data.CellSignal.path;\n');
fprintf(fid,'if nargin<1 || isempty(trange), trange=1:data.CellSignal.t; end\n');
fprintf(fid,'fig = GenFigure_MIPsig_ci(data.CellSignal,trange);\nend\n');
end
