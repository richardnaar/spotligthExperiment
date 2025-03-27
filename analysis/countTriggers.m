function triggerCountsTable = countTriggers(EEG, config)

    % Extract the categories from the config
    categories = fieldnames(config);

    % Create a table to store combinations and their counts
    % Initialize columns based on categories and add a Count column
    vars = [categories; 'Count'];
    combiTable = table('Size', [0 length(vars)], 'VariableTypes', repmat({'cell'}, 1, length(vars)));
    combiTable.Properties.VariableNames = vars;
    
    for eventi = 1:length(EEG.event)
        newRow = table();
        for c = 1:length(categories)
            catName = categories{c};
            fieldVal = EEG.event(eventi).(catName);

            newRow.(catName) = {fieldVal}; % Insert the value as a cell
        end

        newRow.Count = {1}; % Initialize count as 1 in a cell
        
        % Check if a row with the same combination exists
        idx = ismember(combiTable(:, categories), newRow(:, categories));
        
        if any(idx) % if combination exists, increment its count
            countIdx = find(idx);
            combiTable.Count{countIdx} = combiTable.Count{countIdx} + 1;
        else % if combination doesn't exist, add a new row
            combiTable = [combiTable; newRow];
        end
    end

    % Convert the 'Count' column back to double for easier operations
    combiTable.Count = cell2mat(combiTable.Count);
    
    % Return the table
    triggerCountsTable = combiTable;
end


% function triggerCountsTable = countTriggers(EEG, config)
% 
%     % Extract the categories from the config
%     categories = fieldnames(config);
% 
%     % Create a table to store combinations and their counts
%     combiTable = table();
%     combiTable.Combination = {};
%     combiTable.Count = zeros(0,1);
% 
%     for eventi = 1:length(EEG.event)
%         combi = ''; % Initialize a combination string for this event
%         for c = 1:length(categories)
%             catName = categories{c};
%             fieldVal = EEG.event(eventi).(catName);
% 
%             % Build the combination string for this event
%             combi = [combi, catName, ':', fieldVal, ';'];
%         end
% 
%         % Check if combination exists in the table
%         idx = find(strcmp(combiTable.Combination, combi));
%         if isempty(idx) % if combination doesn't exist, add it
%             combiTable = [combiTable; table({combi}, 1, 'VariableNames', {'Combination', 'Count'})];
%         else % if combination exists, increment its count
%             combiTable.Count(idx) = combiTable.Count(idx) + 1;
%         end
%     end
% 
%     % Return the table
%     triggerCountsTable = combiTable;
% end
% 


% function triggerCounts = countTriggers(EEG, config)
% 
%     % Extract the categories from the config
%     categories = fieldnames(config);
% 
%     % Initialize a cell array to store unique combinations of events
%     uniqueCombinations = {};
% 
%     % Go through each event in the EEG structure
%     for eventi = 1:length(EEG.event)
%         combi = ''; % Initialize a combination string for this event
%         for c = 1:length(categories)
%             catName = categories{c};
%             fieldVal = EEG.event(eventi).(catName);
% 
%             % Build the combination string for this event
%             combi = [combi, catName, ':', fieldVal, ';'];
%         end
% 
%         % Add the combination to the uniqueCombinations if it's not already there
%         if ~any(strcmp(uniqueCombinations, combi))
%             uniqueCombinations = [uniqueCombinations; combi];
%         end
%     end
% 
%     % Initialize the counts structure for each unique combination
%     for u = 1:length(uniqueCombinations)
%         sanitizedCombi = strrep(uniqueCombinations{u}, '-', '_'); % replace hyphens with underscores
%         sanitizedCombi = strrep(sanitizedCombi, ':', '_');
%         sanitizedCombi = strrep(sanitizedCombi, ';', '_');
%         triggerCounts.(sanitizedCombi) = 0;
%     end
% 
%     % Go through each event again to count occurrences of each combination
%     for eventi = 1:length(EEG.event)
%         combi = ''; % Initialize a combination string for this event
%         for c = 1:length(categories)
%             catName = categories{c};
%             fieldVal = EEG.event(eventi).(catName);
% 
%             % Build the combination string for this event
%             combi = [combi, catName, ':', fieldVal, ';'];
%         end
% 
%         % Increment count for this combination
%         sanitizedCombi = strrep(combi, '-', '_'); 
%         sanitizedCombi = strrep(sanitizedCombi, ':', '_');
%         sanitizedCombi = strrep(sanitizedCombi, ';', '_');
%         triggerCounts.(sanitizedCombi) = triggerCounts.(sanitizedCombi) + 1;
%     end
% end



% function triggerCounts = countTriggers(EEG, config)
% 
%     % Initialize an empty structure
%     triggerCounts = struct();
% 
%     % Extract the categories from the config
%     categories = fieldnames(config);
% 
%     % Initialize categories in triggerCounts based on config
%     for c = 1:length(categories)
%         catName = categories{c};
%         for v = 1:size(config.(catName), 1)
%             valueName = config.(catName){v, 1};
%             sanitizedValueName = strrep(valueName, '-', '_');  % replace hyphens with underscores
% 
%             % Using dynamic referencing but within known field names from config
%             triggerCounts.(catName).(sanitizedValueName) = 0;
%         end
%     end
% 
%     % Go through each event in the EEG structure
%     for eventi = 1:length(EEG.event)
%         for c = 1:length(categories)
%             catName = categories{c};
%             fieldVal = EEG.event(eventi).(catName);
%             sanitizedFieldVal = strrep(fieldVal, '-', '_');  % again, replace hyphens
% 
%             % Increment count for this category-value pair (knowing it exists from config)
%             triggerCounts.(catName).(sanitizedFieldVal) = triggerCounts.(catName).(sanitizedFieldVal) + 1;
%         end
%     end
% end
