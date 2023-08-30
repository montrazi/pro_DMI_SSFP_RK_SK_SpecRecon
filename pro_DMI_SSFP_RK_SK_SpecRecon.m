% --------------------------------------------------------------------------------------------
%
%     Demo software for RK-SpecRecon and SK-SpecRecon
%             Release version - August 2023
%
% --------------------------------------------------------------------------------------------
%
% The software implements the RK-SpecRecon and SK-SpecRecon algorithm described in the papers:
%
% Montrazi ET, Bao Q, Martinho RP, Peters DC, Harris T, Sasson K, Agemy L, Scherz A, Frydman L. 
% Deuterium imaging of the Warburg effect at sub-millimolar concentrations by joint processing of the kinetic and spectral dimensions.
% NMR Biomed. 2023 Jul 4:e4995. doi: 10.1002/nbm.4995. Epub ahead of print. PMID: 37401393.
%
% --------------------------------------------------------------------------------------------
%
% authors:               Elton Tadeu Montrazi
%                        Qingjia Bao
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

close all; clear; clc;

% %%% ME-SSFP demo %%%
% pathdata=[pwd '\me_ssfp_data_mouse_tumor_pancreas\'];
% expAnatomic=60; %anatomic image folder, # or 'no'
% b0map=7; %fieldmap folder, # or 'no'
% expseq=[8 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59]; %ssfp folders
% timeexp=[0 6.38 13.87 21.37 28.87 36.37 43.87 51.37 58.85 66.33 73.83 81.33 88.85 96.35 103.87...
%     111.38 118.9 126.4 133.92 141.43 148.95 156.45 163.98 171.5 179.03 186.52]; %time for each ssfp (needed for SK-SpecRecon)
% sequence='me-ssfp'; %csi-ssfp or me-ssfp
% %%%%%%%%%%%%%%%%%%%%

%%% CSI-SSFP demo %%%
pathdata=[pwd '\csi_ssfp_data_mouse_tumor_pancreas\'];
expAnatomic=5; %anatomic image folder, # or 'no'
b0map='no'; %fieldmap folder, # or 'no'
expseq=[8 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61]; %ssfp folders
sequence='csi-ssfp'; %csi-ssfp or me-ssfp
%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%% Fitting parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
method='ideal'; % ideal, rk-ideal (RK-SpecRecon), or sk-ideal (SK-SpecRecon)
Nzf=2; %zero filling

%%% RK-SpecRecon %%%
alpha=10; FirstPointReg=[1 0 1]; %zero to remove regularization from the first to second metabolic map 
%%%%%%%%%%%%%%%%%%%%

%%% SK-SpecRecon %%%
if strcmp(method,'sk-ideal')
    K=4; Phi_matrix=Phi_generator_SK_SpecRecon(timeexp,K);
end
%%%%%%%%%%%%%%%%%%%%

imDataParams.PrecessionIsClockwise=-1;
imDataParams.FieldStrength=2.333;
GyromagneticRatio = 42.576;

ppmoff=2.0;
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.7 - ppmoff; % frequency (in ppm)
algoParams.species(1).relAmps = 1; % relative amplitudes
algoParams.species(2).name = 'glucose';
algoParams.species(2).frequency = 3.6 - ppmoff; % frequencies in ppm
algoParams.species(2).relAmps = 1; % relative amplitudes
algoParams.species(3).name = 'lactate';
algoParams.species(3).frequency = 1.2 - ppmoff; % frequencies in ppm
algoParams.species(3).relAmps = 1; % relative amplitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%% Anatomic image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(strcmp(expAnatomic,'no'))
data=OpenBrukerImage([pathdata num2str(expAnatomic) '\pdata\1']);
anaimage=data; Nsl=size(anaimage,3);
anaimage=min(anaimage,max(anaimage,[],'all')/1);
figure; graph1=tiledlayout(3,floor(Nsl/3)+1,'TileSpacing','Compact'); for p=1:Nsl; nexttile; imagesc(anaimage(:,:,p).'); colormap(gray); axis tight; axis off; end
title(graph1,'Anatomic Images');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%% B0 map by 1H %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(strcmp(b0map,'no'))
data=OpenBrukerImage([pathdata num2str(b0map) '\pdata\2']); %Amp FieldMap
AmpB0bruker=data; AmpB0bruker=squeeze(AmpB0bruker(:,:,:,1)); Nsl=size(AmpB0bruker,3); %first echo
AmpB0bruker=permute(AmpB0bruker,[2 1 3]);
% figure; plotAmpfieldmaps=tiledlayout(3,floor(Nsl/3)+1);  title(plotAmpfieldmaps,'2H Field Map from 1H (Amplitude)')
% for p=1:Nsl; nexttile; imagesc(AmpB0bruker(:,:,p)); colormap(jet); caxis([min(AmpB0bruker,[],'all') max(AmpB0bruker,[],'all')]); axis tight; axis off; end; cb=colorbar; cb.Layout.Tile='east';
AmpB0bruker=squeeze(mean(AmpB0bruker,3)); MaskB0=ones(size(AmpB0bruker)); MaskB0(AmpB0bruker<0.3*max(AmpB0bruker,[],'all'))=0;
% figure; imagesc(MaskB0.*AmpB0bruker); colorbar; colormap(jet);

data=OpenBrukerImage([pathdata num2str(b0map) '\pdata\1']); %Read FieldMap
B0bruker=-data/6.5144; Nsl=size(B0bruker,3); %6.5144 -> 1H to 2H
B0bruker=permute(B0bruker,[2 1 3]);
% figure; plotfieldmaps=tiledlayout(3,floor(Nsl/3)+1);  title(plotfieldmaps,'2H Field Map from 1H (Hz)')
% for p=1:Nsl; nexttile; imagesc(B0bruker(:,:,p)); colormap(jet); caxis([min(B0bruker,[],'all') max(B0bruker,[],'all')]); axis tight; axis off; end; cb=colorbar; cb.Layout.Tile='east';

B0bruker=mean(B0bruker(:,:,1:end),3); %Average FieldMap slab
% B0bruker=(B0bruker(:,:,5)); %Central slab

B01H=imboxfilt(MaskB0.*B0bruker); %smooth the FieldMap
figure; imagesc(B01H); colorbar; colormap(jet); title('Average or Central slab - 2H Field Map from 1H (Hz)'); %caxis([-40 40]);
else; B01H=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%% Read CSI-SSFP data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(sequence,'csi-ssfp')
rawdata=OpenBrukerCSI_SSFP([pathdata num2str(expseq(1))]);

Np=rawdata.Np; %Np number of points FID
Nx=rawdata.Nx; %Nx number of phase-encodings
Ny=rawdata.Ny; %Ny number of phase-encodings
tau=1/rawdata.SW; %dwell time
remIni=4;

RowImages=zeros(Nzf*Nx,Nzf*Ny,length(expseq)+1,1,Np-remIni); 
RowKspace=zeros(Nzf*Nx,Nzf*Ny,length(expseq)+1,1,Np-remIni); 
TEs=(remIni+1:Np)*tau;
for jj=1:length(expseq)+1

if jj==length(expseq)+1 %average all data
    RowKspace(:,:,jj,1,:)=sum(RowKspace(:,:,1:end,1,:),3);
    RowImages(:,:,jj,1,:)=sum(RowImages(:,:,1:end,1,:),3);
else

rawdata=OpenBrukerCSI_SSFP([pathdata num2str(expseq(jj))]);
kspace=rawdata.kspace;

kspace=permute(kspace(1+remIni:Np,:,:),[2 3 1]);

RowKspace(Nzf*Nx/2+(-Nx/2+1:Nx/2),Nzf*Ny/2+(-Ny/2+1:Ny/2),jj,1,:)=kspace;
aux=squeeze(RowKspace(:,:,jj,1,:));
aux=fftshift(fft(fft(fftshift(aux),[],1),[],2));
RowImages(:,:,jj,1,:)=aux;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%% Read ME-SSFP data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(sequence,'me-ssfp')
rawdata=OpenBrukerME_SSFP([pathdata num2str(expseq(1))]);

Nx=rawdata.Np; %Np number of points FID
Ny=rawdata.Ny; %Ny number of phase-encodings
TEs=rawdata.TEs; %TEs values
NTEs=length(TEs);

tabexp=expseq;
RowImages=zeros(Nzf*Nx,Nzf*Ny,length(tabexp)+1,1,NTEs); 
RowKspace=zeros(Nzf*Nx,Nzf*Ny,length(tabexp)+1,1,NTEs); 

for nexp=1:length(tabexp)+1

if nexp==length(tabexp)+1 %average all data
    RowKspace(:,:,nexp,1,:)=sum(RowKspace(:,:,1:end,1,:),3);
    RowImages(:,:,nexp,1,:)=sum(RowImages(:,:,1:end,1,:),3);
else
rawdata=OpenBrukerME_SSFP([pathdata num2str(expseq(nexp))]);

kspaceTE=zeros(Nzf*Nx,Nzf*Ny,NTEs); %Make the k-space for each TE
kspaceTE(Nzf*Nx/2+(-Nx/2+1:Nx/2),Nzf*Ny/2+(-Ny/2+1:Ny/2),:)=rawdata.kspace;

imageTE=zeros(Nzf*Nx,Nzf*Ny,NTEs); %Reconstruction images for each TE
for te=1:NTEs; imageTE(:,:,te)=fftn(fftshift(squeeze(kspaceTE(:,:,te)))); end

RowImages(:,:,nexp,1,:)=imageTE;
RowKspace(:,:,nexp,1,:)=kspaceTE;
end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%% Imagem B0 correction for RowImagesTEs %%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:size(RowImages,5); for p=1:size(RowImages,3); RowImages(:,:,p,1,n)=RowImages(:,:,p,1,n).*exp(-1i*pi*2*TEs(n)*B01H); end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%% Fitting IDEAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDataParams.TE=TEs; %TEs (s)
imDataParams.echo_polarity=1;
Nspec=length(algoParams.species);

imDataParams.images=RowImages(:,:,1:end-1,1,:);
forceB0='off';
if strcmp(method,'ideal')
    alpha=0; outParams=RK_SpecRecon_montrazi(imDataParams,algoParams,forceB0,alpha,FirstPointReg);
elseif strcmp(method,'rk-ideal')
    outParams=RK_SpecRecon_montrazi(imDataParams,algoParams,forceB0,alpha,FirstPointReg);
elseif strcmp(method,'sk-ideal')
    outParams=SK_SpecRecon_montrazi(imDataParams,algoParams,forceB0,Phi_matrix);
end

for pp=1:Nspec; AuxSepImage(:,:,:,pp)=outParams.species(pp).amps; end
AuxSepImage(:,:,:,Nspec+1)=outParams.B0map;
SepImage=AuxSepImage; clear AuxSepImage outParams;

%%% Average all data %%%
imDataParams.images=RowImages(:,:,end,1,:);
alpha=0; outParams=RK_SpecRecon_montrazi(imDataParams,algoParams,forceB0,alpha,FirstPointReg);
for pp=1:Nspec; AuxSepImage(:,:,:,pp)=outParams.species(pp).amps; end
AuxSepImage(:,:,:,Nspec+1)=outParams.B0map;
SepImageAve=AuxSepImage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%% Plot SepImages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AuxSepImage=reshape(SepImage,size(SepImage,1),size(SepImage,2)*size(SepImage,3),size(SepImage,4));
if size(SepImage,3)>=15
    AuxSepImage2=AuxSepImage;
    Nparts=round(size(SepImage,3)/2); NpartsF=fix(size(SepImage,3)/2);
    AuxSepImage=zeros(2*size(SepImage,1),Nparts*size(SepImage,2),size(SepImage,4));
    AuxSepImage(1:size(SepImage,1),:,:)=AuxSepImage2(:,1:Nparts*size(SepImage,2),:);
    AuxSepImage(size(SepImage,1)+1:end,1:NpartsF*size(SepImage,2),:)=AuxSepImage2(:,Nparts*size(SepImage,2)+1:end,:);
end
figure; timegraph=tiledlayout(Nspec+1,1,'TileSpacing','Compact');
for p=1:Nspec
    ih=nexttile; imshow(squeeze(abs(AuxSepImage(:,:,p)))); colormap(ih,jet); colorbar;
    [minV,maxV]=bounds(abs(AuxSepImage(:,:,p)),'all'); caxis([minV,maxV]); ylabel(algoParams.species(p).name)
end
ih=nexttile; imshow(AuxSepImage(:,:,end)); colormap(ih,jet); colorbar; [minV,maxV]=bounds(AuxSepImage(:,:,end),'all');
caxis([minV,maxV]); ylabel('B0 map')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%% Plot SepImages Average %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AuxSepImage=SepImageAve;
figure; timegraph=tiledlayout(Nspec+1,1,'TileSpacing','Compact'); title(timegraph,'Average');
for p=1:Nspec
    ih=nexttile; imshow(abs(AuxSepImage(:,:,p))); colormap(ih,jet); colorbar; [minV,maxV]=bounds(abs(AuxSepImage(:,:,p)),'all'); caxis([minV,maxV]);
end
ih=nexttile; imshow(AuxSepImage(:,:,end)); colormap(ih,jet); colorbar; [minV,maxV]=bounds(AuxSepImage(:,:,end),'all'); caxis([minV,maxV]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

