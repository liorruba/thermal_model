function [elementRangeMax, elementRangeMin] = calculateIterationSize(iterationIndex, jobSize, lastJobSize, numberOfJobs, mapSize)
    % It is far easier to work with a zero based index here:
    zeroBasedIndex = iterationIndex - 1;	

    if (zeroBasedIndex == numberOfJobs - 1)
    	elementRangeMax = mapSize^2;
    	elementRangeMin = elementRangeMax - lastJobSize + 1;

    else
    	elementRangeMin = (zeroBasedIndex * jobSize) + 1;
    	elementRangeMax = (zeroBasedIndex + 1) * jobSize;
    end
end