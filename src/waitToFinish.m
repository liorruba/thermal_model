function[] = waitToFinish(secondsToWait, simDir)
	load([simDir, '/config/bin/settings.mat']);

	% If the log file doesn't exist, create it:
	if (~exist([simDir,'/lock']))
		fid=fopen([simDir,'/lock'],'a');
		fclose(fid);
	end

	% Initialize the lock file with a zero.
	lock = 0;
	buff = 1;

	while ( lock ~= settings.numberOfJobs )
		if (buff == 60)
			error('Waiting for longer than 30 minutes, ending program...');
			exit;
		end

		fid=fopen([simDir,'/.lock'],'r');

		if (fid == -1)
			writeToLog(['Cannot open file: ', simDir,'/lock']);
		end

		buff=fscanf(fid,'%d');
		fclose(fid);

		lock = sum(buff);
		writeToLog(['Last job has not finished yet, waiting ', num2str(secondsToWait),' seconds.']);
		pause(secondsToWait);

		buff = buff + 1;
	end
	
	% Remove lock when done:
	if (exist([simDir,'/lock']))
		delete([simDir,'/lock']);
	end
end