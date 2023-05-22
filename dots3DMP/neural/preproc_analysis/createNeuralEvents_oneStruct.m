function S = createNeuralEvents_oneStruct(data)

% pass in data field from one row of datastruct (i.e. one set)

paradigms = fieldnames(data);

allflds = fieldnames(data.(paradigms{1}).pldaps);
numPLDAPSfield = length(allflds);
for par = 1:length(paradigms)

    % kluge
    if ~strcmp(paradigms{par},'RFmapping')
        data.(paradigms{par}).events.headingTheta = nan(size(data.(paradigms{par}).events.trStart));
        data.(paradigms{par}).events.headingPhi   = nan(size(data.(paradigms{par}).events.trStart));
        data.(paradigms{par}).events.hdgOrder     = nan(size(data.(paradigms{par}).events.trStart));
    end

    newflds = fieldnames(data.(paradigms{par}).events);
    allflds = [allflds; setdiff(newflds,allflds)];
end

for fd = 1:length(allflds)
    S.(allflds{fd}) = [];
end

for par = 1:length(paradigms)

    for fd = 1:length(allflds)

        if fd<=numPLDAPSfield
            thisF = data.(paradigms{par}).pldaps.(allflds{fd});

            if strcmp(allflds{fd},'unique_trial_number')
                thisF = num2cell(thisF',[1,size(thisF,1)]);
            end
        else
            try
                thisF = data.(paradigms{par}).events.(allflds{fd});
            catch
                thisF = nan(size(thisF));
            end

        end

        S.(allflds{fd}) = cat(2,S.(allflds{fd}),thisF);

    end

end


