function dydt = cr3bpse(t,y)
% non-dimentional circular restricted 3-body problem for Sun-Earth.  
% y(1) = xp, y(2) = yp, y(3) = zp,
% y(4) = xv, y(5) = yv, y(6) = zv
    xp = y(1);
    yp = y(2);
    zp = y(3);
    xv = y(4);
    yv = y(5);
    zv = y(6);
    
    muSE = 3.036e-6;
    x1 = -muSE;
    x2 = 1-muSE;
    
    xa = xp + 2*yv -(1-muSE)*(xp-x1)/((xp-x1)^2 + yp^2 + zp^2)^(3/2) -muSE*(xp-x2)/((xp-x2)^2 + yp^2 + zp^2)^(3/2);
    ya = yp - 2*xv -(1-muSE)*yp/((xp-x1)^2 + yp^2 + zp^2)^(3/2) -muSE*yp/((xp-x2)^2 + yp^2 + zp^2)^(3/2);
    za = -(1-muSE)*zp/((xp-x1)^2 + yp^2 + zp^2)^(3/2) -muSE*zp/((xp-x2)^2 + yp^2 + zp^2)^(3/2);
    
    dydt = [xv; yv; zv; xa; ya; za];
end