% Calculates and plots the latency between network events and TTL signals.
% Generates an excel file containing event information (event type, timestamp,
% and message).
% Execution time is approx. 1 minute per 10,000 trials, tested on Intel Core m3-6Y30
% Alan Ly, 2018

clear all; clc;

% Input parameters
nTrials = 10000;
saveDirectory = 'C:\OpenEphys\Experiment2';
nFolders = 20;
save2xls = false; % Output results to excel, adds overhead to execution time 

% Preallocate memory for arrays
networkEvents = zeros(2*nTrials*nFolders,1); % seconds
TTL = zeros(2*nTrials*nFolders,1); % seconds
photodiode = zeros(nTrials*nFolders,1); % seconds
trialDuration = zeros(1,1);
intertrialInterval = zeros(1,1); 
pathArray = cell(nFolders, 1); % Cell array used to store save paths

% Initalise structure to hold information regarding 'errors'
e = struct('Description', [], 'Subfolder', [], 'TrialNum', [], 'TTL', [], 'Photodiode', []); 

% Constant variables
SAMPLE_RATE = 30000; % Sample rate of the acquisition board (30kHz)

% Counters
networkEventNo = 1;
TTL_eventNo = 1;
photodiode_eventNo = 1; 
errorNo = 1;

% Determine the path of each sub-folder
dirInfo = dir(saveDirectory);
for i = 1:nFolders 
    folderId = num2str(i);     
    for j = 1:size(dirInfo,1) 
        findFolder = strfind(dirInfo(j).name, [folderId '_2018']);
        if (~isempty(findFolder))
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
    [eventChannel, timestamps, info] = load_open_ephys_data('all_channels.events'); 
    results = [info.eventType timestamps info.eventId]; 

    % For each timestamp in 'all_channels.events,' cross-check with timestamps in 
    % 'messages.events' to find a match. If a match is found, add the message at 
    % that timestamp to corresponding row in the cell array 'msgs.' 
    for i = 1:length(timestamps)
        if (results(i,1) == 5) % Only check network events and not TTL's
            for j = 1:length(parsedTimestamps)
                if (results(i,2) == parsedTimestamps(j))
                    msgs{i,1} = parsedMsgs{j};
                    %break here?
                    break
                end 
            end
         else 
            msgs{i,1} = [];       
        end    
    end 
    
    % Select timestamps of network events corresponding to trial start and complete
    j = 1; 
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
            end
            if (~isempty(findTrialComplete))
                networkEvents(networkEventNo) = results(j,2);
                networkEventNo = networkEventNo + 1; 
                trialCompleteFound = 1;
            end
            j = j+1;
            if (j > length(msgs)) % Fail-safe breaks out of while loop if trial not found
                trialStartFound = 1;
                trialCompleteFound = 1; 
                % Error info here e.g. numTrial missing
                e(errorNo).Description = 'Trial not found'; 
                e(errorNo).Subfolder = pathArray{z};
                e(errorNo).TrialNum = trialId;
                errorNo = errorNo + 1; 
            end 
        end      
    end 

    % Sort digital pulses by channel number
%     networkEventFlag = 1; 
    for i = 1:size(results,1)

        % Separate Daq TTL's (channel #1) into individual array
        if (results(i,1) == 3 && eventChannel(i) == 1) % TTL's eventType = 3
            
            TTL(TTL_eventNo) = results(i,2);
            
            if ((mod(TTL_eventNo, 2) == 1) && results(i,3) ~= 1) %Check for rising edge on odd TTL's               
                e(errorNo).Description = 'TTL mismatch (rising)';
                e(errorNo).Subfolder = pathArray{z}; 
                e(errorNo).TTL = TTL_eventNo; 
                errorNo = errorNo + 1;
            elseif  ((mod(TTL_eventNo, 2) == 0) && results(i,3) ~= 0) %Check for falling edge on even TTL's
                e(errorNo).Description = 'TTL mismatch (falling)';
                e(errorNo).Subfolder = pathArray{z}; 
                e(errorNo).TTL = TTL_eventNo; 
                errorNo = errorNo + 1;
            end 
            
            if (results(i,3) == 1)
                TTL_rising_flag = 1;
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

        % Separate photodiode digital pulses (channel #2) into individal array
%         if (results(i,1) == 5) 
%             networkEventFlag = 1; 
%         end 
%         % If network event has been encountered, log the next photodiode digital pulse 
%         if (results(i,1) == 3 && eventChannel(i) == 2 && networkEventFlag) 
%             photodiode(photodiode_eventNo) = results(i,2);
%             photodiode_eventNo = photodiode_eventNo + 1;
%             networkEventFlag = 0; %Do not log photodiode pulses until another network event is found
%         end
    end 

% Method 2 of finding photodiode pulses
%     TTL_num = find(eventChannel == 1); % Array containing indices of TTL events
%     diff_TTL_num = diff(TTL_num); % Number of array elements between TTL events
%     for i = 1:length(TTL_num) 
%         for j = 1:diff_TTL_num(i)-1
%             if (eventChannel(TTL_num(i+1)) == 2) 
%                 %Log timestamp of photodiode pulse if it follows a TTL
%                 photodiode(photodiode_eventNo) = results(i+1, 2); 
%                 photodiode_eventNo = photodiode_eventNo + 1; 
%             end            
%         end 
%     end 
    
    %Conditions for executing code that generates outputs
    finalFolder = z == length(pathArray); 
    sameLength = length(networkEvents) == length(TTL);
    errorless = isempty(e(1).Description); 

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
        for k = 1:2:length(diff_TTL) 
                trialDuration(1,end+1) = diff_TTL(k,1);
        end 
        
        %Pick out intertrial intevals (every even index)
        for j = 2:2:length(diff_TTL)
            if (mod(j,2*nTrials) ~= 0) %Skip/omit the inter-trial interval between experiments
                intertrialInterval(1,end+1) = diff_TTL(j,1);
            end
        end 
        
        %Remove first zero generated by initialisation
        trialDuration = trialDuration(1,2:end); 
        intertrialInterval = intertrialInterval(1,2:end);
        
        figure(2)
        plot(1:length(trialDuration), trialDuration, 1:length(intertrialInterval), intertrialInterval); 
        title('Trial durations and inter-trial intervals'); 
        xlabel('Trial or inter-trial interval number')
        ylabel('Duration (s)')
        legend('Trial duration', 'Inter-trial interval');
        
        if (isequal(TTL,sort(TTL)))
            disp('Success: TTL''s arrived in order'); 
        else 
            disp('Warning: TTL''s did not arrive in order'); 
        end 
        if (isequal(photodiode, sort(photodiode)))
            disp('Success: Photodiode pulses arrived in order'); 
        else 
            disp('Warning: Photodiode pulses did not arrive in order');
        end 
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