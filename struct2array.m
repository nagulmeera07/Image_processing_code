function array = struct2array(s)
    % Convert structure array to cell array
    if isempty(s)
        array = [];
        return;
    end

    fields = fieldnames(s);
    numFields = numel(fields);
    numStructs = numel(s);

    % Preallocate cell array
    array = cell(numStructs, numFields);

    % Populate cell array
    for i = 1:numFields
        for j = 1:numStructs
            array{j,i} = s(j).(fields{i});
        end
    end
end
