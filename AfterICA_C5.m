%% LEMON: rename the first S210 of each eyes-closed block to "Close"
% This script does not run pop_epoch; it only generates new datasets carrying three "Close" markers.

% Path configuration
inputRoot = '/Users/Re-Re/Desktop/Researches/Alzheimer/MultiCenterVersion/EEG-LEMON/Preprocessed/Data';

if isempty(which('eeglab.m'))
    error('EEGLAB has not been added to the MATLAB path. Please run eeglab first.');
end

% Locate the preprocessed (ICA) datasets
fileList = dir(fullfile(inputRoot, 'sub-*_ICA.set'));

if isempty(fileList)
    error('No sub-*_ICA.set files were found in the input directory.');
end

[~, sortIndex] = sort({fileList.name});
fileList = fileList(sortIndex);

fprintf('Found %d preprocessed datasets.\n', numel(fileList));

% Processing status tracking
participantID = strings(numel(fileList), 1);
status = strings(numel(fileList), 1);
message = strings(numel(fileList), 1);
numberOfCloseMarkers = nan(numel(fileList), 1);
outputFile = strings(numel(fileList), 1);
screenSize = get(0, 'ScreenSize'); 
%% Process each participant
for s = 20:numel(fileList)
    inputFile = fileList(s).name;
    [~, subjectID] = fileparts(inputFile);
    subjectID = erase(subjectID, '_ICA');
    participantID(s) = string(subjectID);
    
    fprintf('[%d/%d] Processing %s.\n', ...
        s, numel(fileList), subjectID);
    
    EEG = pop_loadset( ...
        'filename', inputFile, ...
        'filepath', inputRoot);
    EEG=pop_chanedit(EEG, {'lookup','standard_1005.elc'});
    pop_topoplot(EEG, 0, [1:20] ,fileList(s).name(1:end-4),[4 5] ,0,'electrodes','off');
    Answer0 = inputdlg('Suspect ICs (no more than 4)', 'Suspect', [1 50], {''});
    Answer1 = str2num(Answer0{1});
    close all
    for n = 1:length(Answer1)
        pop_prop( EEG, 0, Answer1(n), NaN, {'freqrange',[2 50] });
        hFig = gcf;
        figWidth = floor((screenSize(3)-10)/4);
        figHeight = 450;
        newPosition = [(n-1)*(figWidth+1), screenSize(4) - figHeight, figWidth, figHeight];
        set(hFig, 'Units', 'pixels', 'Position', newPosition);
    end
    Answer2 = inputdlg('Input ICs to delete', 'Remove ICs', [1 50], {Answer0{1}});
    EEG = pop_subcomp( EEG, str2num(Answer2{1}), 0);
    close all
    EEG.setname=[fileList(s).name(1:end-4),'_Remove_ICs'];
    
    % Normalize event names (used for matching only)
    nEvents = numel(EEG.event);
    eventTypeNormalized = strings(nEvents, 1);
    
    for eventIndex = 1:nEvents
        eventTypeNormalized(eventIndex) = upper(regexprep( ...
            strtrim(string(EEG.event(eventIndex).type)), ...
            '\s+', ''));
    end
    
    % Find all S200 and S210 condition markers
    isClosed = eventTypeNormalized == "S210";
    isOpen = eventTypeNormalized == "S200";
    conditionEventIndices = find(isClosed | isOpen);
    
    if isempty(conditionEventIndices)
        error('%s: no S200 or S210 marker was found.', subjectID);
    end
    
    % Inspect only the condition-marker sequence; ignore boundary, S1, and other events
    conditionTypes = eventTypeNormalized(conditionEventIndices);
    
    % The first S210 of each run of consecutive S210s marks the onset of an eyes-closed block
    isClosedCondition = conditionTypes == "S210";
    isClosedBlockStart = isClosedCondition & ...
        [true; ~isClosedCondition(1:end-1)];
    
    closeEventIndices = ...
        conditionEventIndices(isClosedBlockStart);
    
    % Exactly three eyes-closed blocks must be detected
    if numel(closeEventIndices) ~= 3
        error('%s: detected %d eyes-closed blocks, expected 3.', ...
            subjectID, numel(closeEventIndices));
    end
    
    % Modify only the first S210 of each of the three eyes-closed blocks
    for blockIndex = 1:3
        EEG.event(closeEventIndices(blockIndex)).type = 'Close';
    end
    
    EEG = eeg_checkset(EEG, 'eventconsistency');
    
    % Re-check the number of "Close" markers after modification
    finalEventTypes = strings(numel(EEG.event), 1);
    for eventIndex = 1:numel(EEG.event)
        finalEventTypes(eventIndex) = upper(regexprep( ...
            strtrim(string(EEG.event(eventIndex).type)), ...
            '\s+', ''));
    end
    
    nClose = sum(finalEventTypes == "CLOSE");
    if nClose ~= 3
        error('%s: the number of "Close" markers after modification is %d, expected 3.', ...
            subjectID, nClose);
    end

    fprintf('[%d/%d] %s done: three "Close" markers created.\n', ...
        s, numel(fileList), subjectID);  
    EEG = pop_interp(EEG, [12  16], 'spherical');
    EEG = pop_epoch( EEG, {  'Close'  }, [0  60], 'newname',...
        [fileList(s).name(1:end-4),'_epochs'], 'epochinfo', 'yes');
    EEG = pop_select( EEG, 'channel',{'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'});
    pop_saveset( EEG, 'filename',[fileList(s).name(1:end-4),'_rmvICs.set'],'filepath',...
        '/Users/Re-Re/Desktop/Researches/Alzheimer/MultiCenterVersion/EEG-LEMON/Preprocessed/rmICA');
end
