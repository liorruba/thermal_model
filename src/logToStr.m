% Converts logical variables into string (1 for "true", 0 for "false"):
function [logInStr] = logToStr(logVar)
	if (logVar)
		logInStr = 'true';
	elseif (~logVar)
		logInStr = 'false';
	else
		error('Input is not a logical variable (1 or 0).')
	end
end
