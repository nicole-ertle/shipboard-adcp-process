function [rms_err] = adcp_tide_fit(Pars,Mx_freq,Tstar,Q)
%Calculate RMS error of fit between data of velocity through an
%ADCP underway depth bin and M2 tidal curve of a particular phase and amplitude
%
%Pars: column vector with row 1=residual velocity, followed by alternating pairs of phase(in radians) and amplitude (in meters) of each tidal
%       harmonic in remaining rows (e.g., M2 phase in row 2, M2 amplitude in row 3)
%Mx_freq: vector of frequency of each tidal harmonic, in units of Tstar^-1
%Tstar: fraction of the tidal period from one low tide to next (ranges from 0 to 1)
%Q: bin-averaged velocity vs. Tstar in m/s
%
%Use fminsearch to minimize this function (e.g., adcp_underway_en475_HI.m)
%
%modified from adcp_en475_tide_fit.m
%22 Oct 2018  APH


% %correct size of Pars
% Pars=reshape(Pars,length(Mx_freq)+1,2);

%double length of tidal cycle
Tstar=[Tstar; Tstar+1];
Q=[Q; Q];

%velocity predicted from tidal harmonics 
Q_pred=zeros(size(Q));
for m=1:length(Mx_freq)
    qp=Pars(2*m+1,1).*sin(2*Mx_freq(m)*pi.*Tstar-Pars(2*m,1));
    Q_pred=Q_pred+qp;
end

Q_pred=Q_pred+Pars(1,1); %add residual flux

%RMS error between predicted and observed fluxes
rms_err=sqrt((1/length(Q))*sum((Q-Q_pred).^2));