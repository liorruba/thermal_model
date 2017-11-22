function[] = incrementLock(simDir)
	writeToLog('Incrementing lock.');
	
	fid = fopen([simDir,'/.lock'], 'a');

	% Wait for other processes to stop using the fie:
	counter=1;
	while ( fid < 0 )
		pause(3);
		writeToLog(['Attempting to open lock file once more: ', simDir,'/.lock']);
		fid = fopen([simDir,'/.lock'], 'a');
		counter = counter + 1;
		
		if (counter == 100)
			writeToLog(['Cannot open lock file. Ending simulation. Error code ', fid]);
			error('Simulation ended. Check log for details.');
		end
	end

	% Print the number 1 with a space as a delimiter afterwards.
	fprintf(fid,'%d%s',1,' ');
	fclose(fid);
end

