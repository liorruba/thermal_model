function [jobSize, lastJobSize] = calculateJobSize(mapSize, numberOfJobs)

	% numberOfJobs=2 is a special case, since mod(n,2)=0 or 1.
	if (numberOfJobs == 2)
		jobSize = floor((mapSize^2)./2);
		lastJobSize= ceil((mapSize^2)./2);
	
	else
		jobSize = fix((mapSize^2)./(numberOfJobs - 1));
		lastJobSize= mod(mapSize^2, numberOfJobs - 1);
	end
end