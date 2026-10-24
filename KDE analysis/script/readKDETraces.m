function traces = readKDETraces(cfg)
%READKDETRACES Read one recording/condition per tab, time/value column pairs.
% Observation order, interior missing values and native per-cell axes survive.
% Only the empty padding after an individual column pair is excluded.
% Trace/animal metadata are reconstructed at runtime, without source-file
% registries in the public workbook. SourceColumn refers to the public layout.
assert(isstruct(cfg) && isfield(cfg,'traceLayout'),'kde:InputConfiguration', ...
    'Pass the figure configuration, including traceLayout, to readKDETraces.');
file=cfg.inputFile; sheets=string(sheetnames(file));
[~,base,ext]=fileparts(file); sourceFile=string(base)+string(ext);
required={'TraceID','Condition','FlyIndex','CellIndex','CellType','CellName', ...
    'SourceFile','SourceSheet','SourceColumn','SourceFormat','IsReference'};
traces=repmat(struct('meta',table(),'time',[],'signal',[]),0,1);
if cfg.traceLayout=="aligned_pairs"
    assert(numel(sheets)==numel(cfg.conditions),'kde:ConditionCount','Expected one tab per configured ZT condition.');
elseif cfg.traceLayout=="fly_cells"
    conditionCounts=zeros(size(cfg.conditions));
    starts=nan(size(cfg.conditions));
    for c=1:numel(starts)
        token=regexp(char(cfg.conditions(c)),'ZT(\d+)','tokens','once');
        assert(~isempty(token),'kde:ConditionName','A condition needs its nominal ZT.');
        starts(c)=str2double(token{1});
    end
elseif cfg.traceLayout=="continuous_cells"
    validateattributes(cfg.legacySheetCount,{'numeric'},{'scalar','integer','nonnegative','<=',numel(sheets)});
end
for s=1:numel(sheets)
    raw=readcell(file,'Sheet',s);
    assert(size(raw,1)>=2 && mod(size(raw,2),2)==0,'kde:RecordingShape', ...
        'A recording tab must contain time/value column pairs: %s',sheets(s));
    header=string(raw(1,:)); names=erase(header(1:2:end),'_Time_h');
    assert(all(endsWith(header(1:2:end),'_Time_h')) && ...
        all(header(2:2:end)==names+"_Value") && numel(unique(names))==numel(names), ...
        'kde:RecordingHeaders','Expected <cell>_Time_h / <cell>_Value headers.');
    values=nan(size(raw,1)-1,size(raw,2)); cells=raw(2:end,:);
    for k=1:numel(cells)
        value=cells{k};
        if isnumeric(value) && isscalar(value),values(k)=double(value);
        else
            assert(isempty(value) || all(ismissing(value),'all'),'kde:RecordingNumeric', ...
                'Recording observations must be numeric or empty.');
        end
    end
    isE=startsWith(lower(names),'5th') | ~cellfun(@isempty,regexp(cellstr(names),'^E\d+$','once'));
    condition="";fly=0;format="";reference=0;
    switch cfg.traceLayout
        case "aligned_pairs"
            condition=cfg.conditions(s);format="paired";
        case "native_single"
            token=regexp(char(sheets(s)),'CT(\d+)','tokens','once');
            assert(~isempty(token) && numel(names)==1,'kde:NativeTab', ...
                'Native-trace tabs must retain a CT label and one time/value pair.');
            condition="CT"+string(token{1});format="single";
        case "fly_cells"
            nominal=nan(size(names));
            for j=1:numel(names)
                t=values(:,2*j-1);first=find(isfinite(t),1);
                assert(~isempty(first),'kde:MissingTime','A cell has no valid time values.');
                nominal(j)=round(t(first));
            end
            assert(all(nominal==nominal(1)),'kde:ConditionTime','Cell start times imply different ZT groups.');
            c=find(starts==nominal(1));
            assert(isscalar(c),'kde:ConditionTime','Recording start ZT does not match a configured condition.');
            conditionCounts(c)=conditionCounts(c)+1;
            condition=cfg.conditions(c);fly=conditionCounts(c);format="cell";
        case "continuous_cells"
            condition="continuous";fly=s;
            if s<=cfg.legacySheetCount,format="old";else,format="new";end
            assert(sum(isE)<=1,'kde:CellOrder','A recording can contain at most one E cell.');
            reference=find(isE,1);
            if isempty(reference),reference=1;end
        otherwise
            error('kde:InputLayout','Unknown trace layout: %s',cfg.traceLayout);
    end
    for j=1:numel(names)
        pair=values(:,2*j+[-1 0]);
        last=find(any(~isnan(pair),2),1,'last');
        assert(~isempty(last),'kde:EmptyTrace','A column pair has no observations.');
        pair=pair(1:last,:);
        signal=pair(:,2);finite=signal(isfinite(signal));
        assert(all(finite==0 | finite==1),'kde:NonBinaryInput','The binary workbook contains nonbinary signal values.');
        switch cfg.traceLayout
            case "aligned_pairs",type="unspecified";
            case "fly_cells",type=cfg.cellType;
            otherwise
                if isE(j),type="E";else,type="M";end
        end
        index=numel(traces)+1;
        traceID=string(sprintf('T%03d',index));
        metadata=table(traceID,string(condition),double(fly),double(j),string(type),names(j), ...
            sourceFile,sheets(s),double(2*j),string(format),double(j==reference), ...
            'VariableNames',required);
        traces(index,1)=struct('meta',metadata,'time',pair(:,1),'signal',pair(:,2));
    end
end
end
