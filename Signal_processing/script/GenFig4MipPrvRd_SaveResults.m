function GenFig4MipPrvRd_SaveResults(~,fig)
% Save to the recording currently being reviewed, including after relocation.
disp('Saving...');
i=fig.UserData.CellIndex;
root=fullfile(fig.UserData.path,fig.UserData.SigInfo.signalfolder);
[~,~,files]=listSignalFolders(root,i);
CellSignal=fig.UserData;
save(files{1},'CellSignal','-append');
disp(['cell_' num2str(i,'%02d') '_tracked_signal.mat saved']);
end
