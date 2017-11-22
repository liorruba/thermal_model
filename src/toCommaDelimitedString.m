function[string] = toCommaDelimitedString(vector)
% This function converts a vector into a comma delimited string.
%
% Input variables:
% vector - the vector about to be converted.
%
% Output variables:
% string - the comma seperated string.
%
% Written by Lior Rubanenko, Weizmann Institute of Science, Jan 2015
%

string = sprintf('%.02f,' , vector);
% Remove the last comma:
string = string(1:end-1);
