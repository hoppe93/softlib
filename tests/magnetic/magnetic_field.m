function [B,Babs,gradB,curlB] = magnetic_field(x0,B0,Rm,zm,sigmaB,sigmaI,q,rDqDr)

    x = x0(1);
    y = x0(2);
    z = x0(3);
    
    R = hypot(x, y);
    r = hypot(z-zm, Rm-R);
    
    sintheta = (z-zm) ./ r;
    costheta = (Rm-R) ./ r;
    sinphi = y ./ R;
    cosphi = x ./ R;
    
    rhat     = [-costheta*cosphi;-costheta*sinphi;sintheta];
    thetahat = [sintheta*cosphi; sintheta*sinphi; costheta];
    phihat   = [-sinphi; cosphi; 0];
    
    B = B0*Rm ./ R * (sigmaI*r/(q*Rm)*thetahat + sigmaB*phihat);
    Babs = norm(B);
    
    DBr = B0*Rm./R * (r./(q*Rm).^2 .* (1-rDqDr)./sqrt((r/(q*Rm)).^2 + 1) + costheta./R.*sqrt((r/(q*Rm)).^2+1));
    DB0 =-B0*Rm*sintheta / R.^2 .* sqrt((r/(q*Rm))^2+1);
    
    gradB = DBr * rhat + DB0 * thetahat;
    
    curlB = sigmaI * B0./(q.*R) .* (Rm./R + 1 - rDqDr./q) * phihat;

end
