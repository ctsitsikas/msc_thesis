%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: PDF.m
%Author: Christos Tsitsikas
%Email: ctsitsikas@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = PDF(x, loc, scale) % takes as input x values, location and scale parameter, and gives Gumbel PDF

for i =1:length(x)
    y(i) = (1./scale).*exp((((loc-x(i))./scale)-exp((loc-x(i))./scale)));
end
end