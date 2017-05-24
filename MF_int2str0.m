function str = MF_int2str0(int, dig)
% _
% Convert integer to string of specified length
% FORMAT str = MF_int2str0(int, dig)
% 
%     int - integer
%     dig - digits
% 
%     str - string
% 
% FORMAT str = MF_int2str0(int, dig) converts the integer int to a string
% using int2str.m and then adds zeros at the left until the resulting
% string has a length of dig. Example: MF_int2str0(42,5) = '00042'.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/03/2015, 07:55 (V0.3/V10)
%  Last edit: 26/03/2015, 09:20 (V0.3/V10)


% Convert to string
%-------------------------------------------------------------------------%
str = int2str(int);

% Add zeros at the left
%-------------------------------------------------------------------------%
while length(str) < dig
    str = strcat('0',str);
end;