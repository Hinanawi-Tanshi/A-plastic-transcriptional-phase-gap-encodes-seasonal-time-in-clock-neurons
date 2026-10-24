function NormalizeResults_spline250822_all_cell_min_20x20(SigInfo,i)
% Historical filename retained: baseline correction uses rloess (span 0.2),
% not a spline. Original formulas below are preserved.
%% input parameters
tvalue = SigInfo.tvalue;
zvalue = SigInfo.zvalue;
xvalue = SigInfo.xvalue;
yvalue = SigInfo.yvalue;
signalsize = SigInfo.signalsize;
bgsize = SigInfo.bgsize;
signalfolder = SigInfo.signalfolder;
p = SigInfo.p;
[~,~,files] = listSignalFolders(fullfile(p,signalfolder),i);
filename = files{1};
data = load(filename,'CellSignal');
CellSignal = data.CellSignal;
assert(size(CellSignal.MIPsig,2)==7 && ~isfield(CellSignal,'NormalizedMIPsig'), ...
    'MS2:NormalizedInput','Normalization already exists; preserve it and start from pre-normalization results.');
%%
[~,T]=size(CellSignal.MIP);
MIPsig = CellSignal.MIPsig;
CellSignal.MIPsigBeforeNormalization = MIPsig; % retain reviewed measurements
sigflag=MIPsig(:, end);
normalizedMIPsig = sigflag;

bg_smooth = smooth(1:T, MIPsig(:, 3), 0.2, 'rloess');

if size(MIPsig, 2) ~= 8
    MIPsig(:, end+1) = MIPsig(:, end);
end

all_signal = 1;

for t = 1:T
    if sigflag(t) == 0
        all_signal = 0;
        MIPsig(t, end-1) = MIPsig(t, 3);
        MIPsig(t, 1) = MIPsig(t, 3) - bg_smooth(t) + MIPsig(1, 3);   % bgi / mean_cell
        if MIPsig(t, 1) < 0
            MIPsig(t, 1) = 0;
        end
    else 
        MIPsig(t, end-1) = MIPsig(t, 2);
        MIPsig(t, 1) = MIPsig(t, 2) - bg_smooth(t) + MIPsig(1, 3);   % signali / mean_cell
    end
end

normalizedMIPsig(:, 2) = MIPsig(:, 1);

normalizedBg = MIPsig(:, 3) - bg_smooth + MIPsig(1, 3);

% if all_signal == 1 || all_signal == 0
%     min_normalized_res = min(normalizedBg);
%     mean_normalized_res = mean(normalizedBg, "all");
% else
%    min_normalized_res = min(normalizedMIPsig(:, 2));
%    mean_normalized_res = mean(normalizedMIPsig(:, 2), "all"); 
% end 

min_normalized_res = min(normalizedMIPsig(:, 2));
mean_normalized_res = mean(normalizedMIPsig(:, 2), "all");
min_normalized_bg = min(normalizedBg);
mean_normalized_bg = mean(normalizedBg, "all");


% raw delta F / F : min & mean from de_back result
normalizedMIPsig(:, 3) = (normalizedMIPsig(:, 2) - min_normalized_res) / mean_normalized_res;

% 20x20 delta F / F : min & mean from de_back bg
normalizedMIPsig(:, 4) = (normalizedMIPsig(:, 2) - min_normalized_bg) / mean_normalized_bg;

% bg delta F / F
normalizedMIPsig(:, 5) = (normalizedBg - min(normalizedBg)) / mean(normalizedBg);

% normalizedMIPsig(:, 4) = MIPsig(:, 2) - bg_smooth;

CellSignal.NormalizedMIPsig = normalizedMIPsig;
CellSignal.MIPsig = MIPsig;
%%
save(filename, 'CellSignal','-append');
disp(['cell ' num2str(CellSignal.CellIndex,'%02d') ' is done.'])