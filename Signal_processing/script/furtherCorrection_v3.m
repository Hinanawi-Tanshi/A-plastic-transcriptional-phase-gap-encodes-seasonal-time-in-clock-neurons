function furtherCorrection_v3(SigInfo,i)
% Recalculate OFF-frame measurements near the mean ON-frame position.
% The existing ON/OFF flags and KDE confidence thresholds are retained.
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
    'MS2:NormalizedInput','OFF-frame correction must precede normalization.');
%%
[~,T]=size(CellSignal.MIP);
MIPsig = CellSignal.MIPsig;
sigflag=MIPsig(:,end);
sig_r = [];sig_c = [];
if sum(sigflag)~=0
    sig_r = round(mean(MIPsig(sigflag==1,5)));
    sig_c = round(mean(MIPsig(sigflag==1,6)));
else
    sz = [];
    for img = 1:5:T
        mipsize = size(CellSignal.MIP{img});
        sz = [sz;mipsize];
    end
    sz = round(mean(sz)/2);
    sig_r = sz(1);
    sig_c = sz(2);
end

for t = 1:T
    if sigflag(t) == 0
        cellMIP = CellSignal.MIP{t};
        [sig,row,col] = GetSignal_v4((cellMIP),signalsize,bgsize,[sig_r sig_c]);    % changed
        MIPsig(t,1:end-1) = [sig,row,col];
    end
end
CellSignal.MIPsig = MIPsig;
%%
save(filename, 'CellSignal','-append');
disp(['cell ' num2str(CellSignal.CellIndex,'%02d') ' is done.'])
% SaveMFile_MIP([cellsignalfolder '/cell_' num2str(CellSignal.CellIndex,'%02d') '_ProofReading.m']);
% SaveMFile_MIPsig([cellsignalfolder '/cell_' num2str(CellSignal.CellIndex,'%02d') '_ProofReading_MIPsig.m']);
%SaveMFile_MIPsig_ci([cellsignalfolder '/cell_' num2str(CellSignal.CellIndex,'%02d') '_ProofReading_MIPsig_ci.m']);