% Calculates and plots the latency of network events 
% Generates an excel file containing event information (event type, timestamp,
% event channel, rising/falling edge and message).
% Execution time is approx. 1 minute per 10,000 trials, tested on Intel Core m3-6Y30
% Alan Ly, 2018

clear all; close all; clc;

% Input parameters
nTrials = 1000;
nFolders = 100;
saveDirectory = 'C:\OpenEphys\Experiment5-long\08082018'; % Parent folder under which Open Ephys recording sub-folders are located
save2xls = false; % Output results to excel, adds overhead to execution time 
importNSData = true; % Import Neurostim data and perform clock sync 
nsFilename = 'results08082018.mat'; % Filename of .mat file containing Neurostim experiment data

% Preallocate memory for arrays
networkEvents = zeros(2*nTrials*nFolders,1); % seconds
TTL = zeros(2*nTrials*nFolders,1); % seconds
photodiode = zeros(nTrials*nFolders,1); % seconds
trialDuration = zeros(nTrials*nFolders,1); % seconds
intertrialInterval = zeros(nTrials*nFolders - nFolders, 1); % seconds
ns_Photodiode = zeros(nTrials*nFolders,1); % ms 
ns_NetworkEvents = zeros(2*nTrials*nFolders,1); % ms 
ns_TTL = zeros(2*nTrials*nFolders,1);  % ms 
pathArray = cell(nFolders, 1); % Cell array used to store save paths

% Initalise structure to hold information regarding 'errors'
e = struct('Description', [], 'Subfolder', [], 'TrialNum', [], 'TTL', [], 'Photodiode', []); 

% Constant variables
SAMPLE_RATE = 30000; % Sample rate of the acquisition board (30kHz)

% Counters
networkEventNo = 1;
TTL_eventNo = 1;
photodiode_eventNo = 1; 
ns_networkEventNo = 1; 
ns_TTL_eventNo = 1;
ns_photodiode_eventNo = 1;
trialDurationNo = 1; 
intertrialIntervalNo = 1; 
errorNo = 1;

%Initialise variables
rising_IndexParity = 1; 
falling_IndexParity = 0;

% Determine the path of each sub-folder
dirInfo = dir(saveDirectory);
for i = 1:nFolders 
    folderId = num2str(i);     
    for j = 1:size(dirInfo,1) 
        findFolder = strfind(dirInfo(j).name, [folderId '_2018']);
        if (findFolder == 1)
            pathArray{i, 1} = [saveDirectory  '\' dirInfo(j).name];            
        end 
    end 
end  

for z = 1:length(pathArray) % Repeat code below for each sub-directory
    clearvars events msgs
    cd(pathArray{z,1}); 

    % Extract event messages and their timestamps from the 'messages.events' file 
    fid = fopen('messages.events'); 

    if (fid == -1) 
        error('Could not open ''messages.events'' file. Please check that the file is on MATLAB''s search path');
    end   

    % Determine the size of the file in bytes     
    fseek(fid, 0, 'eof'); 
    fileSize = ftell(fid); 
    fseek(fid, 0, 'bof');     

    eventNo = 1; 
    eventLength = 1;

    % Keep reading file until the file pointer reaches the end of the file.
    while (ftell(fid) < fileSize)   
        % Read a single byte 
        % char*1 preserves the encoding scheme of the file (8-bit ASCII)
        events(eventNo, eventLength) = fread(fid, 1, 'char*1');    

        % If newline (event) detected, start a new row 
        % '10' in ASCII-8 corresponds to a newline
        if (events(eventNo, eventLength) == 10)     
            eventLength = 1; % Reset column counter 
            eventNo = eventNo + 1; % Increment row counter
        else 
            eventLength = eventLength + 1; 
        end 
    end 
    
    fclose(fid); 
    data = char(events); % Convert 8-bit ASCII to characters. 

    parsedTimestamps = zeros(size(data,1), 1); 
    parsedMsgs = cell(size(data,1), 1); 

    % Separate timestamps and messages into individual arrays
    for i = 1:size(data,1) 
        % The separator between timestamp and message is a space.
        spaceIndex = regexp(data(i,:),' ');  

        % Convert timestamp to seconds by dividing by the sample rate
        parsedTimestamps(i) = str2double(data(i,1:spaceIndex(1)))/SAMPLE_RATE; 

        parsedMsgs{i} = data(i, spaceIndex(1):end);
    end 

    % Output the timestamps of TTL's and network events stored in 'all_channels.events'
    % Function is available at analysis-tools repo of Open Ephys
    [eventChannel, timestamps, info] = load_open_ephys_data('all_channels.events'); 
    results = [info.eventType timestamps info.eventId]; 

    % For each timestamp in 'all_channels.events,' cross-check with timestamps in 
    % 'messages.events' to find a match. If a match is found, add the message at 
    % that timestamp to corresponding row in the cell array 'msgs.' 
    j = 1; 
    jLog = 1; 
    for i = 1:length(timestamps)
        msgFound = 0; 
        if (results(i,1) == 5) % Only check network events and not TTL's
            while (~msgFound) % an assumption is made that messages arrive in order
                if (results(i,2) == parsedTimestamps(j))
                    msgs{i,1} = parsedMsgs{j};
                    msgFound = 1;
                    j = j + 1; 
                    jLog = j;
                else 
                    j = j + 1; 
                end 
                if (j > length(parsedTimestamps)) % If match not found, break out of loop
                    msgFound = 1; 
                    j = jLog + 1; 
                end 
            end
        else 
            msgs{i,1} = [];       
        end    
    end 
    
    % Select timestamps of network events corresponding to trial start and complete
    j = 1; % Counter used to increment through msgs{j}
    jLog = 1; % used to log 'j' when a trial is found  
    for i = 1:nTrials
        trialId = num2str(i);
        trialStartFound = 0;
        trialCompleteFound = 0; 
        while (~(trialStartFound && trialCompleteFound) )
            findTrialStart = strfind(msgs{j}, ['_T' trialId '_']);
            findTrialComplete = strfind(msgs{j}, ['l' trialId 'c']);
            if (~isempty(findTrialStart))
                networkEvents(networkEventNo) = results(j,2);
                networkEventNo = networkEventNo + 1; 
                trialStartFound = 1;
                jLog = j; 
            end
            if (~isempty(findTrialComplete))
                networkEvents(networkEventNo) = results(j,2);
                networkEventNo = networkEventNo + 1; 
                trialCompleteFound = 1;
                jLog = j; 
            end
            j = j+1;
            if (j > length(msgs)) % Fail-safe breaks out of while loop if trial not found
                trialStartFound = 1;
                trialCompleteFound = 1;
                j = jLog; % Reset the msg{j} counter to the last known 'hit'
                % Error info here e.g. numTrial missing
                e(errorNo).Description = 'Trial not found'; 
                e(errorNo).Subfolder = pathArray{z};
                e(errorNo).TrialNum = trialId;
                errorNo = errorNo + 1; 
            end 
        end      
    end 

    % Sort digital pulses by channel number
    TTL_rising_flag = 0; % Initialise TTL_rising_flag to prevent error when TTL is missing 
    for i = 1:size(results,1)

        % Separate Daq TTL's (channel #1) into individual array
        if (results(i,1) == 3 && eventChannel(i) == 1) % TTL's eventType = 3
            
            TTL(TTL_eventNo) = results(i,2);
            
            if ((mod(TTL_eventNo, 2) == rising_IndexParity) && results(i,3) ~= 1) %Check for TTL's with rising edge              
                e(errorNo).Description = 'Missing TTL (rising)';
                e(errorNo).Subfolder = pathArray{z}; 
                e(errorNo).TTL = TTL_eventNo; 
                errorNo = errorNo + 1;
                rising_IndexParity = ~rising_IndexParity; % flip the parity (odd or even) of the index on which rising edges are looked for in TTL
                falling_IndexParity = ~falling_IndexParity; 
            elseif  ((mod(TTL_eventNo, 2) == falling_IndexParity) && results(i,3) ~= 0) %Check for TTL's with falling edge
                e(errorNo).Description = 'Missing TTL (falling)';
                e(errorNo).Subfolder = pathArray{z}; 
                e(errorNo).TTL = TTL_eventNo; 
                errorNo = errorNo + 1;
                rising_IndexParity = ~rising_IndexParity; 
                falling_IndexParity = ~falling_IndexParity;
            end 
            
            if (results(i,3) == 1) % eventId == 1 denotes rising edge
                TTL_rising_flag = TTL_rising_flag + 1;
            end
            
            if (TTL_rising_flag > 1 && results(i,3) == 1)
                e(errorNo).Description = 'Missing photodiode';
                e(errorNo).Subfolder = pathArray{z}; 
                e(errorNo).Photodiode = TTL_eventNo; % beforeFrame TTL at which photodiode went missing
                errorNo = errorNo + 1; 
            end 
            
            TTL_eventNo = TTL_eventNo + 1;

        end
        
        %Separate photodiode digital pulses (channel #2) into individual array.
        %Only take the first pulse immediately following a TTL with rising edge
        if (results(i,1) == 3 && eventChannel(i) == 2 && TTL_rising_flag) 
            
            photodiode(photodiode_eventNo) = results(i,2);
            
            if (results(i,3) ~= 1) %Check that the photodiode pulse is rising edge
                e(errorNo).Description = 'Photodiode mismatch';
                e(errorNo).Subfolder = pathArray{z}; 
                e(errorNo).Photodiode = photodiode_eventNo; 
                errorNo = errorNo + 1;
            end 
            
            TTL_rising_flag = 0; %Do not log photodiode pulses until another TTL is found
            photodiode_eventNo = photodiode_eventNo + 1;
            
        end
    end  

    
    %Conditions for executing code that generates outputs
    finalFolder = z == length(pathArray); 
    sameLength = length(find(networkEvents ~= 0)) == length(find(TTL ~= 0)); 
    errorless = isempty(e(1).Description); 
    outputState = 0; 
    if (finalFolder && sameLength && errorless) %Only compute outputs if the outer loop is on the final folder 

        %Caculate and plot message-TTL latencies
        latency = abs(networkEvents - TTL); 

        figure(1)
        plot(1:length(networkEvents), latency*10^3);
        title('message-TTL latencies')
        xlabel('Event number'); 
        ylabel('Latency (ms)'); 

        %Calculate and plot trialDuration and intertrial intervals
        diff_TTL = diff(TTL); 
        
        %Pick out trialDurations (every odd index) 
        for i = 1:2:length(diff_TTL) 
                trialDuration(trialDurationNo) = diff_TTL(i,1);
                trialDurationNo = trialDurationNo + 1; 
        end 
        
        %Pick out intertrial intevals (every even index)
        for j = 2:2:length(diff_TTL)
            if (mod(j,2*nTrials) ~= 0) %Skip/omit the inter-trial interval between experiments
                intertrialInterval(intertrialIntervalNo) = diff_TTL(j,1);
                intertrialIntervalNo = intertrialIntervalNo + 1; 
            end
        end 
        
        figure(2)
        plot(1:length(trialDuration), trialDuration*10^3, 1:length(intertrialInterval), intertrialInterval*10^3); 
        title('Trial durations and inter-trial intervals'); 
        xlabel('Trial or inter-trial interval number')
        ylabel('Duration (ms)')
        legend('Trial duration', 'Inter-trial interval');                       
        
        if (~isequal(TTL,sort(TTL)))
        disp('Warning: TTL''s did not arrive in order');               
        end
        if (~isequal(photodiode, sort(photodiode)))
            disp('Warning: Photodiode pulses did not arrive in order');                    
        end
        
        if (importNSData) 
            load([saveDirectory '\' nsFilename], 'paramLog', 'latencyStruct'); % Implement exist()  
            
            % Extract neurostim times for stimulus onset
            trialNum = 1; 
            for i = 1:length(paramLog(1).data)
                if (paramLog(1).trialTime(i) > 0 && paramLog(1).trial(i) == trialNum) % photodiode pulses occur after trial start
                   ns_Photodiode(ns_photodiode_eventNo) = paramLog(1).time(i);
                   ns_photodiode_eventNo = ns_photodiode_eventNo + 1; 
                   trialNum = trialNum + 1;
                   if (trialNum == nTrials + 1)
                       trialNum = 1; 
                   end 
                end 
            end 
            
             if (ns_photodiode_eventNo ~= nTrials*nFolders + 1)
                e(errorNo).Description = 'Missing Neurostim photodiode';                
                e(errorNo).TrialNum = ns_photodiode_eventNo; 
                errorNo = errorNo + 1;
             end 
            
            % Extract neurostim times for network events
            trialNum = 1; 
            signalEdge = 1;
            instancesFound = 0; 
            for i = 1:length(paramLog(2).data)
                if (paramLog(2).data{i} == signalEdge && paramLog(2).trial(i) == trialNum) 
                    ns_NetworkEvents(ns_networkEventNo) = paramLog(2).time(i);
                    ns_networkEventNo = ns_networkEventNo + 1;
                    signalEdge = ~signalEdge;
                    instancesFound = instancesFound + 1; 
                    if (instancesFound == 2) 
                        trialNum = trialNum + 1;
                        instancesFound = 0; 
                    end
                    if (trialNum == nTrials + 1) 
                        trialNum = 1; 
                    end 
                end 
            end 
            
            if (ns_networkEventNo ~= 2*nTrials*nFolders + 1) 
                e(errorNo).Description = 'Missing Neurostim network event';
                e(errorNo).TrialNum = ns_networkEventNo; 
                errorNo = errorNo + 1; 
            end 
            
            % Extract neurostim times for TTL's
            trialNum = 1; 
            signalEdge = 1;
            instancesFound = 0; 
            for i = 1:length(paramLog(3).data)
                if (paramLog(3).data{i} == signalEdge && paramLog(3).trial(i) == trialNum) % Implement error checking...whether rise, fall, rise ,fall
                    ns_TTL(ns_TTL_eventNo) = paramLog(3).time(i);
                    ns_TTL_eventNo = ns_TTL_eventNo + 1; 
                    signalEdge = ~signalEdge;
                    instancesFound = instancesFound + 1; 
                    if (instancesFound == 2) 
                        trialNum = trialNum + 1;
                        instancesFound = 0; 
                    end 
                    if (trialNum == nTrials + 1) 
                        trialNum = 1; 
                    end 
                end 
            end 
            
           if (ns_TTL_eventNo ~= 2*nTrials*nFolders + 1) 
               e(errorNo).Description = 'Missing Neurostim TTL'; 
               e(errorNo).TrialNum = ns_TTL_eventNo; 
               errorNo = errorNo + 1; 
           end 
           
            % Re-align events within each trial to the time of the photodiode pulse
            offsetArray = repelem(photodiode,2);
            networkEvents_offset = networkEvents - offsetArray;
            TTL_offset = TTL - offsetArray;     
                       
            ns_offsetArray = repelem(ns_Photodiode,2);
            ns_networkEvents_offset = (ns_NetworkEvents - ns_offsetArray)*10^-3; 
            ns_TTL_offset = (ns_TTL - ns_offsetArray)*10^-3; 
 
            
            figure(3)
            ns_latency = abs(networkEvents_offset - ns_networkEvents_offset)*10^3;
            plot(ns_latency); 
            title('Latency of network messages'); 
            xlabel('Event no.'); 
            ylabel('Latency (ms)')
            
            
            figure(4) 
            plot(abs(TTL_offset - ns_networkEvents_offset)*10^3);
            title('Time difference between Neurostim network event and Open Ephys TTL')
            xlabel('Event no.'); 
            ylabel('Latency (ms)');
            
            figure(5)
            plot(abs(TTL_offset - ns_TTL_offset)*10^3); 
            title('TTL latency');
            xlabel('Event no.'); 
            ylabel('Latency (ms)'); 
            
            figure(6)
            plot(abs(((ns_Photodiode - ns_TTL(2:2:end))*10^-3) - (photodiode - TTL(2:2:end))))
            title('Photodiode latency'); 
            xlabel('Event no.'); 
            ylabel('Latency (s)');
        end 
        
        outputState = 1; 
    end 
    
    if (~outputState && finalFolder) % give user feedback in the command window 
        disp('Errors encountered...plots will not be generated');
    elseif (~isempty(e(1).Description) && finalFolder)
        disp('Error encountered');
    elseif (finalFolder) 
        disp('Analysis completed...'); 
    end 
    
    if (save2xls) 
        if (exist('events.xlsx', 'file') ~= 0) 
            delete events.xlsx
        end 

        %Write the results array to an Excel file for improved readability. 
        filename = 'events.xlsx';
        colHeaders = {'EventType', 'eventChannel', 'Rising/Falling', 'Timestamps (s)', 'Message'}; 
        xlswrite(filename, colHeaders, 'Sheet1', 'A1');
        xlswrite(filename, results(:,1), 'Sheet1', 'A2');
        xlswrite(filename, eventChannel, 'Sheet1', 'B2'); 
        xlswrite(filename, results(:,3), 'Sheet1' , 'C2');
        xlswrite(filename, results(:,2), 'Sheet1', 'D2');
        xlswrite(filename, msgs, 'Sheet1', 'E2');
    end 
end 