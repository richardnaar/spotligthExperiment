function EEG = gtrigf(EEG, config, visualCheck, newTrialIndicator)

    % Populate the 'trig' structure based on the config structure
    categories = fieldnames(config);
    trig = struct();
    
    for c = 1:length(categories)
        catName = categories{c};
        values = config.(catName);
        
        % Create category entries
        trig.(catName) = values(:, 1:2);
        
        % Assign pin information directly from values
        trig.(['pin_' catName]) = values(:, 3);
    end

    trialn = -1;  % Start counting from zero

    % Loop through all the events
    for eventi = 1:length(EEG.event)
        if EEG.event(eventi).type ~= 1024  
            trigger = dec2bin(EEG.event(eventi).type);  % Convert to binary

    % Process each category
    catNames = fieldnames(trig);
    for i = 1:length(catNames)
        if ~contains(catNames{i}, 'pin_')
            pinRange = trig.(['pin_' catNames{i}]){1};
            fieldValue = trigger(pinRange(1):pinRange(end));
            if ~any(strcmp(fieldValue, trig.(catNames{i})(:,2)))
                warning(['Unrecognized trigger value: ', fieldValue, ' for category: ', catNames{i}]);
                EEG.event(eventi).(catNames{i}) = 'UNKNOWN';
            else
                EEG.event(eventi).(catNames{i}) = trig.(catNames{i}){ strcmp(fieldValue, trig.(catNames{i})(:,2)), 1};
            end
        end
    end


            % Advance trial counter if it is the start of the trial 
            if strcmp(EEG.event(eventi).(catNames{1}), newTrialIndicator)
                trialn = trialn + 1;
                EEG.event(eventi).('thisN') = trialn;
            end

            EEG.event(eventi).thisN = trialn;

            if visualCheck
                eventDisplay = [];
                for c = 1:length(categories)
                    eventDisplay = [eventDisplay, EEG.event(eventi).(categories{c}), '-'];
                end
                fprintf([eventDisplay, '\n']);
                fprintf([trigger, '\n']);
            end

        else
            for c = 1:length(categories)
                EEG.event(eventi).(categories{c}) = 'START';
            end
        end
    end
end
