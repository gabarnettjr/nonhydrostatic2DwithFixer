function [ cellAreas, areaRatios ] = getCellAreas( xx, zz )

x1 = xx(1:end-1,1:end-1);  z1 = zz(1:end-1,1:end-1);
x2 = xx(2:end,1:end-1);    z2 = zz(2:end,1:end-1);
x3 = xx(2:end,2:end);      z3 = zz(2:end,2:end);
x4 = xx(1:end-1,2:end);    z4 = zz(1:end-1,2:end);

cellAreas = x1.*z2-z1.*x2 + x2.*z3-z2.*x3 + x3.*z4-z3.*x4 + x4.*z1-z4.*x1;
cellAreas = abs(cellAreas(:)) ./ 2;

areaRatios = cellAreas ./ sum(cellAreas);