% Generate the laser somewhere else, the LGS reflects it to the telescope

% Let's give LUVOIR a laser that is 10X bigger and more powerful than LGS
pwr_laser_scope = 10*pwr_laser;
D_laser_scope = 10*D_laser;

LGS_face_area = 0.2*0.3; % 20x30 cm face for the retro
eff_retro_d = sqrt(4*LGS_face_area/pi);

% LUVOIR -> LGS
[PrxTx,PhotrxTx,appMagTx,bwTx] = linkbudgetG(pwr_laser_scope,D_laser_scope,range_LGS,lambda,eff_retro_d);

% LGS back to LUVOIR
[PrxRx,PhotrxRx,appMagRx,bwRx] = linkbudgetG(PrxTx,eff_retro_d,range_LGS,lambda,scope_d);

% Case 2: using LUVOIR's main telescope itself as the transmitter!

% LUVOIR -> LGS
[PrxTx2,PhotrxTx2,appMagTx2,bwTx2] = linkbudgetG(pwr_laser_scope,scope_d,range_LGS,lambda,eff_retro_d);

% LGS back to LUVOIR
[PrxRx2,PhotrxRx2,appMagRx2,bwRx2] = linkbudgetG(PrxTx2,eff_retro_d,range_LGS,lambda,scope_d);

% Giant ABL-TMT facility (Is such a thing even possible? Yes it is.)

pwr_laser_gnd = 1e6;
D_laser_gnd = 30;
lambda_abl = 1315e-9;
range_L2 = ((muSE/3)^(1/3))*AU;

% ABL-TMT -> LGS, assuming no atmosphere (!!!)
[PrxTx3,PhotrxTx3,appMagTx3,bwTx3] = linkbudgetG(pwr_laser_gnd,D_laser_gnd,range_L2,lambda_abl,eff_retro_d);

% LGS back to LUVOIR
[PrxRx3,PhotrxRx3,appMagRx3,bwRx3] = linkbudgetG(PrxTx3,eff_retro_d,range_LGS,lambda_abl,scope_d);

% Now let's take atmospheric effects into account
D_laser_eff_gnd = 0.1; % Fried parameter ~10 cm eff D

% ABL-TMT -> LGS
[PrxTx4,PhotrxTx4,appMagTx4,bwTx4] = linkbudgetG(pwr_laser_gnd,D_laser_eff_gnd,range_L2,lambda_abl,eff_retro_d);

% LGS back to LUVOIR
[PrxRx4,PhotrxRx4,appMagRx4,bwRx4] = linkbudgetG(PrxTx4,eff_retro_d,range_LGS,lambda_abl,scope_d);