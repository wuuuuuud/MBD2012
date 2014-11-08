function [ output ] = getTimeString( )
%getTimeString return an escaped datetime string
output=strrep(datestr(clock),':','-');

end

