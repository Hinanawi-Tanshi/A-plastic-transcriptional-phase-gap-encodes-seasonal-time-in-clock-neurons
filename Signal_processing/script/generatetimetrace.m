function timetrace = generatetimetrace (CellSignal)
%% generate time trace
Firstframe = datevec(CellSignal.T0); % ZT time for first captured frame
FrameRate = CellSignal.FrameInterval ;% in second
%==========================end of input parameter==========================
FrameInterval = ((FrameRate)/60)/60; % in hour
timevec_hzt = Firstframe(4);
timevec_m = Firstframe(5);
timevec_num = timevec_hzt + timevec_m/60;
timetrace = timevec_num: FrameInterval: (timevec_num + FrameInterval * (CellSignal.t-1));
% validTlength = size(CellSignal.MIP,2);
% timetrace = timetrace(1:validTlength);
end