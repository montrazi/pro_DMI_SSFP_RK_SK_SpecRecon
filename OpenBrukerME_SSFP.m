function rawdata=OpenBrukerME_SSFP(pathfile)

acqpBruker=fopen([pathfile '\acqp']); p=1; infoBruker=[]; infoBruker{1}=fgetl(acqpBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(acqpBruker); end
fclose(acqpBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1)
    if contains(infoBruker{p},'##$GO_raw_data_format='); VisuCoreWordType=erase(infoBruker{p},'##$GO_raw_data_format='); end
    if contains(infoBruker{p},'##$ACQ_O1_list=')
        ACQ_O1_list=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; ACQ_O1_list=[ACQ_O1_list str2num(infoBruker{pp})]; pp=pp+1; end
    end
end
switch VisuCoreWordType
    case 'GO_32BIT_SGN_INT'; precision='int32';
    case 'GO_16BIT_SGN_INT'; precision='int16';
    case 'GO_32BIT_FLOAT'; precision='float32';
end

recoBruker=fopen([pathfile '\pdata\1\reco']); p=1; infoBruker=[]; infoBruker{1}=fgetl(recoBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(recoBruker); end
fclose(recoBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1)
    if contains(infoBruker{p},'##$RECO_rotate='); dimRot=erase(infoBruker{p},'##$RECO_rotate='); dimRot=erase(dimRot,'('); dimRot=str2num(erase(dimRot,')')); end
    if contains(infoBruker{p},'##$RECO_rotate=')
        RECO_rotate=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; RECO_rotate=[RECO_rotate str2num(infoBruker{pp})]; pp=pp+1; end
    end
end
RECO_rotate=reshape(RECO_rotate,dimRot(1),dimRot(2));
RECO_rotate=RECO_rotate(:,1);

methodBruker=fopen([pathfile '\method']); p=1; infoBruker=[]; infoBruker{1}=fgetl(methodBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(methodBruker); end
fclose(methodBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1)
    if contains(infoBruker{p},'##$DelayArrayForTE=')
        TEs=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; TEs=[TEs str2num(infoBruker{pp})]; pp=pp+1; end
    end
    if contains(infoBruker{p},'##$PVM_Matrix='); AcqMatrix=[str2num(infoBruker{p+1}) str2num(infoBruker{p+2})]; end
    if contains(infoBruker{p},'##$PVM_NRepetitions='); NR=str2num(erase(infoBruker{p},'##$PVM_NRepetitions=')); end
    if contains(infoBruker{p},'##$PVM_DigNp='); Np=str2num(erase(infoBruker{p},'##$PVM_DigNp=')); end %Np number of points FID
    if contains(infoBruker{p},'##$PVM_DigNp='); Np=str2num(erase(infoBruker{p},'##$PVM_DigNp=')); end
end
TEs=TEs/1000;

pathBruker=fopen([pathfile '\fid']); data=fread(pathBruker,precision); fclose(pathBruker); %read fid

data=data(1:2:end)+1i*data(2:2:end); %complex data
Nbuf=128+fix(Np/128)*128;
data=reshape(data,Nbuf,[]); data=data(1:Np,:);
data=reshape(data,Np,length(TEs)*AcqMatrix(1),NR);

if NR~=1; data=squeeze(sum(data,3));end

NTEs=length(TEs);
kspace=zeros(Np,AcqMatrix(2),NTEs);
for te=1:NTEs; kspace(:,:,te)=squeeze(data(:,te:NTEs:end)); end

for te=1:NTEs
auxdata=squeeze(kspace(:,:,te));
auxdata=conj(auxdata);

dims=[size(auxdata,1),size(auxdata,2),size(auxdata,3),size(auxdata,4)];
phase_matrix=ones(size(auxdata));
for index=1:length(size(auxdata))
    f=1:dims(index);
    phase_vektor=exp(1i*2*pi*RECO_rotate(index)*f);
    switch index
       case 1
           phase_matrix=phase_matrix.*repmat(phase_vektor', [1, dims(2), dims(3), dims(4)]);
       case 2
           phase_matrix=phase_matrix.*repmat(phase_vektor, [dims(1), 1, dims(3), dims(4)]);
       case 3
           tmp(1,1,:)=phase_vektor;
           phase_matrix=phase_matrix.*repmat(tmp, [dims(1), dims(2), 1, dims(4)]);
       case 4
           tmp(1,1,1,:)=phase_vektor;
           phase_matrix=phase_matrix.*repmat(tmp, [dims(1), dims(2), dims(3), 1]);
    end   
    clear phase_vekor f tmp
end
kspace(:,:,te)=auxdata.*phase_matrix;
end

for te=1:NTEs; kspace(:,:,te)=kspace(:,:,te)*exp(-1i*pi*2*TEs(te)*ACQ_O1_list); end %Slice phase correction

rawdata.Nx=AcqMatrix(2);
rawdata.Ny=AcqMatrix(1);
rawdata.Np=Np;
rawdata.TEs=TEs;
rawdata.kspace=kspace;

