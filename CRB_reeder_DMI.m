close all; clear; clc;

% --------------------------------------------------------------------------------------------
%
%              Optimization for DMI using ME-SSFP and CSI-SSFP
%           (TE increment/dwell time choice, maximization of NSA)
%                    Release version - August 2023
%
% --------------------------------------------------------------------------------------------
%
% Optimization for DMI using ME-SSFP and CSI-SSFP when B0 is known:
%
% Peters DC, Markovic S, Bao Q, Preise D, Sasson K, Agemy L, Scherz A, Frydman L.
% Improving deuterium metabolic imaging (DMI) signal-to-noise ratio by spectroscopic multi-echo bSSFP: A pancreatic cancer investigation.
% Magn Reson Med. 2021 Nov;86(5):2604-2617. doi: 10.1002/mrm.28906. Epub 2021 Jun 30. PMID: 34196041.
%
% Reeder SB, Wen Z, Yu H, Pineda AR, Gold GE, Markl M, Pelc NJ. 
% Multicoil Dixon chemical species separation with an iterative least-squares estimation method.
% Magn Reson Med. 2004 Jan;51(1):35-45. doi: 10.1002/mrm.10675. PMID: 14705043.
%
% For a complete discussion with B0 is also unknown:
%
% Pineda AR, Reeder SB, Wen Z, Pelc NJ.
% Cramér-Rao bounds for three-point decomposition of water and fat. 
% Magn Reson Med. 2005 Sep;54(3):625-35. doi: 10.1002/mrm.20623. PMID: 16092102.
%
% More about ME-SSFP sequence:
%
% Leupold J, Månsson S, Petersson JS, Hennig J, Wieben O.
% Fast multiecho balanced SSFP metabolite mapping of (1)H and hyperpolarized (13)C compounds.
% MAGMA. 2009 Aug;22(4):251-6. doi: 10.1007/s10334-009-0169-z. Epub 2009 Apr 15. PMID: 19367422.
%
% --------------------------------------------------------------------------------------------
%
% authors:               Elton Tadeu Montrazi
%                        Lucio Frydman
%
% web page:              https://www.weizmann.ac.il/chembiophys/Frydman_group/software
%                        https://github.com/montrazi                                           
%
% contact:               elton.montrazi@gmail.com
%                        lucio.frydman@weizmann.ac.il
%
% --------------------------------------------------------------------------------------------
% Copyright (c) Weizmann Institute of Science.
% All rights reserved.
% This work should be used for nonprofit purposes only.
% --------------------------------------------------------------------------------------------


B0=15.2; LarFreq=6.5357; %B0(T), Larmor's frequency (MHz/T)
flipPulse1=pi*(60/180); %Flip Angle
TR=0.0115; %TR to plot

ppmoff=2.0; %carrier frequency (ppm)

algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.7 - ppmoff; % frequency (in ppm)
algoParams.species(1).relAmps = 1; % relative amplitudes
algoParams.species(1).T1 = 0.4; % T1 (s)
algoParams.species(1).T2 = 0.04; % T2 (s)

algoParams.species(2).name = 'glucose';
algoParams.species(2).frequency = 3.7 - ppmoff; % frequencies in ppm
algoParams.species(2).relAmps = 1; % relative amplitudes
algoParams.species(2).T1 = 0.053; % T1 (s)
algoParams.species(2).T2 = 0.053; % T2 (s)

algoParams.species(3).name = 'lactate';
algoParams.species(3).frequency = 1.2 - ppmoff; % frequencies in ppm
algoParams.species(3).relAmps = 1; % relative amplitudes
algoParams.species(3).T1 = 0.27; % T1 (s)
algoParams.species(3).T2 = 0.27; % T2 (s)

NTEs=5; %N echoes
t=linspace(1,10000,100)/10^6;


%%%%%%%%%%%%%%%% Cramér–Rao bound %%%%%%%%%%%%%%%%%%
Nspec=length(algoParams.species);
for p=1:Nspec; wfreq(p,1)=B0*LarFreq*algoParams.species(p).frequency; end

A=zeros(NTEs,length(wfreq));
NSA=zeros(length(t),length(wfreq));
for p=1:length(t)
    TEs(:,1)=(1:NTEs)*t(p);
    %%%% Make the matrix A %%%%
    for j=1:NTEs
        A(j,:)=exp(1i*2*pi*wfreq*TEs(j,1))';
    end
    NSA(p,:)=1./diag(inv(A'*A));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; plot(t*1000,NSA);
xlabel('TE increment (ms)'); ylabel('NSA'); legend(algoParams.species.name);

figure; plot(t*1000*NTEs,NSA);
xlabel('Acq. Time (ms)'); ylabel('NSA'); legend(algoParams.species.name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% SSFP signal %%%%%%%%%%%%%%%%%%%%
for p=1:Nspec
    wtr(1,p)=2*pi*B0*LarFreq*algoParams.species(p).frequency;
    T1(1,p)=algoParams.species(p).T1;
    T2(1,p)=algoParams.species(p).T2;
    AmpSa(1,p)=algoParams.species(p).relAmps;
end

%%%%%% Amplitude at TE=TR/2 in function of TR %%%%%%%%%%%%%%
alphapi=1*pi; %pi for alternated flip pulse (alpha,-alpha) and 0 without alternating
figAmpTR=figure; hold on;
ttr=linspace(0,30,1000)/1000;
for p=1:length(wtr)
    E1=exp(-ttr./T1(p)); E2=exp(-ttr./T2(p)); theta=ttr.*wtr(p)+alphapi;
    D=(1-E1.*cos(flipPulse1)).*(1-E2.*cos(theta))-(E1-cos(flipPulse1)).*(E2-cos(theta)).*E2;
    MyTR=(1-E1).*((1-E2.*cos(theta)).*sin(flipPulse1))./D; MxTR=(1-E1).*(E2.*sin(flipPulse1).*sin(theta))./D;
    AmpTR=AmpSa(p)*(MxTR+1i*MyTR).*exp(-1i*theta); [val,idx]=min(abs(ttr-TR));
    plot(ttr,abs(AmpTR),'--',ttr(idx),abs(AmpTR(idx)),'pr'); text(ttr(idx),abs(AmpTR(idx)),num2str(abs(AmpTR(idx))));
    pleng(p)=plot(ttr,abs(AmpTR)); pleng2=plot(ttr(idx),abs(AmpTR(idx)),'pr');
    text(ttr(idx),abs(AmpTR(idx)),num2str(abs(AmpTR(idx))));
end
xlabel('TR (s)'); ylabel('Amplitude'); legend([pleng pleng2],{algoParams.species.name,['TR=' num2str(TR)]});
hold off;

%%%%%% Amplitude in function of alpha %%%%%%%%%%%%%%
figure; hold on; alpha=(pi/180)*linspace(0,180,181);
for p=1:length(wtr)
    E1=exp(-TR./T1(p)); E2=exp(-TR./T2(p)); theta=TR*wtr(p)+alphapi;
    D=(1-E1.*cos(alpha)).*(1-E2.*cos(theta))-(E1-cos(alpha)).*(E2-cos(theta)).*E2;
    MyTR=(1-E1).*((1-E2.*cos(theta)).*sin(alpha))./D; MxTR=(1-E1).*(E2.*sin(alpha).*sin(theta))./D;
    AmpTR=AmpSa(p)*(MxTR+1i*MyTR).*exp(-1i*theta);
    AmpTRCRB(p)=abs(AmpTR((flipPulse1/pi)*180+1));
    plotAmpTR(p,1:2)=plot(alpha*(180/pi),abs(AmpTR),flipPulse1*(180/pi),abs(AmpTR((flipPulse1/pi)*180+1)),'pr');
    pleng(p)=plot(alpha*(180/pi),abs(AmpTR)); pleng2=plot(flipPulse1*(180/pi),abs(AmpTR((flipPulse1/pi)*180+1)),'pr');
    text(flipPulse1*(180/pi),abs(AmpTR((flipPulse1/pi)*180+1)),num2str(abs(AmpTR((flipPulse1/pi)*180+1))));
end
xlabel('alpha (^o)'); ylabel('Amplitude');legend([pleng pleng2],{algoParams.species.name,[num2str(flipPulse1*(180/pi)) '^o']});
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

