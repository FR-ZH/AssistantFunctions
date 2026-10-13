eeglab
close all
clear
clc
Folder = '/Users/Re-Re/Desktop/Researches/Alzheimer/PreProcess/Rest2';   %Read ICA file
SaveFolder = '/Users/Re-Re/Desktop/Researches/Alzheimer/AfterICA/Rest2';     %Save 
setFiles = dir(fullfile(Folder, '*.set'));
screenSize = get(0, 'ScreenSize'); 
for i = 3:length(setFiles)
    workFile = setFiles(i).name;
    EEG = pop_loadset('filename',workFile,'filepath',[Folder,'/']);
    pop_topoplot(EEG, 0, [1:15] ,workFile(1:end-4),[4 5] ,0,'electrodes','on');
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
    EEG.setname=[workFile(1:end-4),'_Remove_ICs'];
    pop_saveset( EEG, 'filename',[workFile(1:end-4),'_rmvICs.set'],'filepath',[SaveFolder,'/']);
    clear EEG Answer0 Answer1 Answer2 n hFig newPosition figWidth figHeight ans workFile
end