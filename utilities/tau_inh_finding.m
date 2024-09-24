function tauMargin = tau_inh_finding(M)

eg = eig(M);
phase = atan2(imag(eg),real(eg));
phase = phase + (phase<0)*2*pi;
om = abs(eg);
id_nonzero = om>1e-10;
tauMargin = min((phase(id_nonzero)-pi/2)./om(id_nonzero));  % t_inh margin
tauMargin = floor(tauMargin*1e4)/1e4;
% tauinhDM = pi/2;
if any((abs(phase(id_nonzero))-pi/2)<0)
    error('Delay-free system is unstable')
end


end