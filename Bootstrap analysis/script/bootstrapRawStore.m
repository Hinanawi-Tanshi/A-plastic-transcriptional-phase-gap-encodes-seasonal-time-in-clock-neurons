function result = bootstrapRawStore(operation,location,argument)
% Fig2C input: binary and continuous-intensity workbooks, one recording per tab.
% Match YYMMDD FLYn tab names across workbooks. Derive nominal start CT from
% the first finite timestamp, rounded to the nearest hour without modulo 24.
% Existing analysis loaders use virtual recording paths supplied by this store.
persistent cache
location=char(location);
stop=regexpi(location,'\.xlsx','end','once');
assert(~isempty(stop),'bootstrap:RawWorkbookPath','An XLSX input path is required.');
binaryPath=location(1:stop);
relative=string(regexprep(strrep(location(stop+1:end),'\','/'),'^/+',''));
if strcmp(operation,'configure')
    cache=configure(binaryPath,argument);
    result=[];
    return
end
assert(~isempty(cache) && strcmp(cache.path,binaryPath), ...
    'bootstrap:RawNotConfigured','Configure the split recording inputs before analysis.');
[bytes,modified]=fileSignature(cache.path,cache.options.intensityFile);
if ~isequal(bytes,cache.bytes) || ~isequal(modified,cache.modified)
    cache=configure(binaryPath,cache.options);
end
switch operation
    case 'isfolder'
        result=any(cache.files.Batch==relative | cache.folderKeys==relative);
    case 'folders'
        names=unique(cache.files.Folder(cache.files.Batch==relative),'sorted');
        result=makeEntries(names,location,true);
    case 'files'
        names=sort(cache.files.FileName(cache.folderKeys==relative));
        result=makeEntries(names,location,false);
    case {'sheets','read'}
        use=cache.fileKeys==relative;
        assert(sum(use)==1,'bootstrap:RawFileKey','Recording is missing or ambiguous.');
        index=cache.files.FileID(use);
        if isempty(cache.values{index})
            [names,data,ct]=readRecording(cache.files.PhysicalFile(index),cache.files.SheetIndex(index));
            assert(ct==cache.files.StartCT(index),'bootstrap:RawTimeMismatch', ...
                'Binary/intensity start CT differs for recording %s.',cache.files.FileName(index));
            cache.cellNames{index}=names;
            cache.values{index}=data;
        end
        if strcmp(operation,'sheets')
            result=cache.cellNames{index};
        else
            cellIndex=find(strcmp(cache.cellNames{index},char(argument)));
            assert(isscalar(cellIndex),'bootstrap:RawCellKey','Cell is missing or ambiguous.');
            data=cache.values{index}(:,2*cellIndex+[-1 0]);
            result=table(data(:,1),data(:,2),'VariableNames',{'Time_h','Value'});
        end
    otherwise
        error('bootstrap:RawOperation','Unknown operation: %s',operation);
end
end

function cache=configure(binaryPath,config)
assert(isstruct(config) && isfield(config,'intensityFile') && isfield(config,'batchSheetCounts'), ...
    'bootstrap:RawConfiguration','Specify intensityFile and batchSheetCounts.');
options=struct('intensityFile',char(config.intensityFile),'batchSheetCounts',config.batchSheetCounts);
[bytes,modified]=fileSignature(binaryPath,options.intensityFile);
binaryNames=string(sheetnames(binaryPath)); binaryNames=binaryNames(:);
intensityNames=string(sheetnames(options.intensityFile)); intensityNames=intensityNames(:);
binaryKeys=recordingKeys(binaryNames);
intensityKeys=recordingKeys(intensityNames);
counts=options.batchSheetCounts(:);
validateattributes(counts,{'numeric'},{'vector','integer','positive','finite'});
assert(sum(counts)==numel(binaryNames),'bootstrap:RawBatchCount', ...
    'batchSheetCounts must account for every binary recording sheet.');
orphan=~ismember(intensityKeys,binaryKeys);
assert(~any(orphan),'bootstrap:OrphanIntensity', ...
    'Intensity tabs have no matching binary recording: %s',strjoin(intensityNames(orphan),', '));
[paired,intensityIndex]=ismember(binaryKeys,intensityKeys);
nBinary=numel(binaryNames); nFiles=nBinary+numel(intensityNames);
FileID=(1:nFiles).'; Batch=strings(nFiles,1); Folder=strings(nFiles,1);
FileName=strings(nFiles,1); PhysicalFile=strings(nFiles,1);
SheetIndex=zeros(nFiles,1); StartCT=zeros(nFiles,1);
cache=struct('path',binaryPath,'options',options,'bytes',bytes,'modified',modified);
cache.values=cell(nFiles,1); cache.cellNames=cell(nFiles,1);
batchEnds=cumsum(counts); row=0;
for i=1:nBinary
    [names,data,ct]=readRecording(binaryPath,i);
    binaryValues=data(:,2:2:end);
    finite=binaryValues(isfinite(binaryValues));
    assert(all(finite==0 | finite==1),'bootstrap:NonBinaryInput', ...
        'Binary workbook contains nonbinary values in tab %s.',binaryNames(i));
    batch=find(i<=batchEnds,1);
    id=string(sprintf('Recording%04d',i));
    % Retain the original recording order and the original paired-only scope
    % of the intensity analysis. Unpaired binary records still enter KDE.
    if paired(i), kinds=[2 1]; else, kinds=3; end
    for kind=kinds
        row=row+1;
        Batch(row)="Batch"+string(batch); Folder(row)="CT"+string(ct);
        StartCT(row)=ct;
        if kind==2
            FileName(row)=id+"_intensity.xlsx";
            PhysicalFile(row)=string(options.intensityFile);
            SheetIndex(row)=intensityIndex(i);
        else
            if kind==1,FileName(row)=id+"_binary.xlsx";
            else,FileName(row)=id+"_unpaired.xlsx";end
            PhysicalFile(row)=string(binaryPath); SheetIndex(row)=i;
            cache.cellNames{row}=names; cache.values{row}=data;
        end
    end
end
assert(row==nFiles,'bootstrap:RawCount','Recording count mismatch.');
cache.files=table(FileID,Batch,Folder,FileName,PhysicalFile,SheetIndex,StartCT);
cache.folderKeys=Batch+"/"+Folder;
cache.fileKeys=cache.folderKeys+"/"+FileName;
assert(numel(unique(cache.fileKeys))==nFiles,'bootstrap:RawFileKey','Duplicate virtual recording key.');
end

function keys=recordingKeys(names)
keys=strings(size(names));
for k=1:numel(names)
    token=regexpi(char(strtrim(names(k))),'^(\d{6})\s+FLY\s*(\d+)$','tokens','once');
    assert(~isempty(token),'bootstrap:RecordingName', ...
        'Use recording tab names in YYMMDD FLYn format: %s',names(k));
    keys(k)=string(token{1})+" fly"+string(str2double(token{2}));
end
assert(numel(unique(keys))==numel(keys),'bootstrap:DuplicateRecordingName', ...
    'A workbook contains duplicate date/FLY recording keys.');
end

function [names,numeric,ct]=readRecording(path,index)
raw=readcell(path,'Sheet',index);
assert(mod(size(raw,2),2)==0 && size(raw,1)>=2,'bootstrap:RawShape', ...
    'Each recording needs time/value column pairs and at least one observation.');
headers=string(raw(1,:));
cellNames=erase(headers(1:2:end),'_Time_h');
assert(all(endsWith(headers(1:2:end),'_Time_h')) && ...
    all(headers(2:2:end)==cellNames+"_Value") && numel(unique(cellNames))==numel(cellNames), ...
    'bootstrap:RawHeaders','Expected <cell>_Time_h / <cell>_Value column pairs.');
numeric=nan(size(raw,1)-1,size(raw,2)); data=raw(2:end,:);
for k=1:numel(data)
    value=data{k};
    if isnumeric(value) && isscalar(value), numeric(k)=double(value);
    else
        assert(isempty(value) || all(ismissing(value),'all'), ...
            'bootstrap:RawNumeric','Recording observations must be numeric or empty.');
    end
end
starts=nan(1,numel(cellNames));
for k=1:numel(cellNames)
    times=numeric(:,2*k-1); first=find(isfinite(times),1);
    assert(~isempty(first),'bootstrap:MissingRecordingTime','A cell has no valid timestamps.');
    starts(k)=round(times(first));
end
assert(all(starts==starts(1)),'bootstrap:RawTimeMismatch', ...
    'Cells in a recording disagree on nominal start CT.');
ct=starts(1); names=cellstr(cellNames(:));
end

function [bytes,modified]=fileSignature(binaryPath,intensityPath)
a=dir(binaryPath); b=dir(intensityPath);
assert(isscalar(a) && isscalar(b),'bootstrap:MissingRawWorkbook', ...
    'Both binary and continuous-intensity workbooks are required.');
bytes=[a.bytes b.bytes]; modified=[a.datenum b.datenum];
end

function entries=makeEntries(names,folder,isDirectory)
entries=repmat(struct('name','','folder',folder,'isdir',isDirectory),numel(names),1);
for k=1:numel(names),entries(k).name=char(names(k));end
end
