function[] = writeToLog(inputText, displayMsg)
	% This function writes the given data to the log file.
	%
	% inputText - The text to be written.
	% displayMsg - True to display entry on screen, false otherwise.
	% simDir - The simulation directory.
	%

	global settings;

	if (nargin == 1)
		displayMsg = false;
	end

	% Open log file to writing:
	logFid = fopen([settings.dirPath.logs,'log.txt'],'a');
    
    if logFid == -1
        error('Cannot create log');
    end

	% Add date and MATLAB "stamp" to string:
	inputText = [datestr(clock), '\t', inputText, '\n'];

	if (displayMsg)
		fprintf(1, inputText);
	end

	fprintf(logFid, inputText);
	if (fclose(logFid) == -1)
        error('Cannot properly close log file');
    end    
    
end