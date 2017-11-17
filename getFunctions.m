function [ phi, phi_x, phi_y, phiHV ] = getFunctions( rbfType )

if strcmp( rbfType, 'mq' )
	phi   = @(rad,x,y,param)      mq(rad,x,y,param);
	phi_x = @(rad,x,y,param)    mq_x(rad,x,y,param);
	phi_y = @(rad,x,y,param)    mq_y(rad,x,y,param);
	phiHV = @(rad,x,y,param,K)  mqHV(rad,x,y,param,K);
elseif strcmp( rbfType, 'iq' )
	phi   = @(rad,x,y,param)      iq(rad,x,y,param);
	phi_x = @(rad,x,y,param)    iq_x(rad,x,y,param);
	phi_y = @(rad,x,y,param)    iq_y(rad,x,y,param);
	phiHV = @(rad,x,y,param,K)  iqHV(rad,x,y,param,K);
elseif strcmp( rbfType, 'phs' )
	phi   = @(rad,x,y,param)      phs(rad,x,y,param);
	phi_x = @(rad,x,y,param)    phs_x(rad,x,y,param);
	phi_y = @(rad,x,y,param)    phs_y(rad,x,y,param);
	phiHV = @(rad,x,y,param,K)  phsHV(rad,x,y,param,K);
elseif strcmp( rbfType, 'ga' )
	phi   = @(rad,x,y,param)      ga(rad,x,y,param);
	phi_x = @(rad,x,y,param)    ga_x(rad,x,y,param);
	phi_y = @(rad,x,y,param)    ga_y(rad,x,y,param);
	phiHV = @(rad,x,y,param,K)  gaHV(rad,x,y,param,K);
end
