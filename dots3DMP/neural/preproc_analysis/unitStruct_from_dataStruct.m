function unitStruct = unitStruct_from_dataStruct(dataStruct, par)

% convert dataStruct, which is organized by recording set, into 'unitStruct', with one row per unit
% this might be more intuitive for running certain analyses

% TODO, add information about unit channel and depth...

unitStruct = struct();

T = 0;
for s = 1:length(dataStruct)

    if ~isfield(dataStruct(s).data,par), continue, end

    try
        units    = dataStruct(s).data.(par).units;
    catch
        continue
    end

    numUnits = length(units.cluster_id);

    fldnames = fieldnames(dataStruct(s));
    fldnames(strcmp(fldnames,'data')) = [];

    ufldnames = fieldnames(units);
    
    for u=1:numUnits
        
        for ff = 1:length(fldnames)
            unitStruct(T+u).(fldnames{ff}) = dataStruct(s).(fldnames{ff});
        end


        for ff = 1:length(ufldnames)
            unitStruct(T+u).(ufldnames{ff}) = units.(ufldnames{ff})(u);
        end

        unitStruct(T+u).events = dataStruct(s).data.(par).events;

    end

    T = length(unitStruct);
end
