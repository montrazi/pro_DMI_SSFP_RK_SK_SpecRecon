function rawdata=OpenBrukerCSI_SSFP(pathfile)

acqpBruker=fopen([pathfile '\acqp']); p=1; infoBruker=[]; infoBruker{1}=fgetl(acqpBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(acqpBruker); end
fclose(acqpBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1); if contains(infoBruker{p},'##$GO_raw_data_format='); VisuCoreWordType=erase(infoBruker{p},'##$GO_raw_data_format='); end; end
switch VisuCoreWordType
    case 'GO_32BIT_SGN_INT'; precision='int32';
    case 'GO_16BIT_SGN_INT'; precision='int16';
    case 'GO_32BIT_FLOAT'; precision='float32';
end

methodBruker=fopen([pathfile '\method']); p=1; infoBruker=[]; infoBruker{1}=fgetl(methodBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(methodBruker); end
fclose(methodBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1)
    if contains(infoBruker{p},'##$AverageList=')
        AverageList=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; AverageList=[AverageList str2num(infoBruker{pp})]; pp=pp+1; end
    end
    if contains(infoBruker{p},'##$PVM_Matrix='); AcqMatrix=[str2num(infoBruker{p+1}) str2num(infoBruker{p+2})]; end
    if contains(infoBruker{p},'##$PVM_NRepetitions='); NR=str2num(erase(infoBruker{p},'##$PVM_NRepetitions=')); end
    if contains(infoBruker{p},'##$PVM_DigSw='); SW=str2num(erase(infoBruker{p},'##$PVM_DigSw=')); end
    if contains(infoBruker{p},'##$PVM_DigNp='); Np=str2num(erase(infoBruker{p},'##$PVM_DigNp=')); end %Np number of points FID
end

pathBruker=fopen([pathfile '\rawdata.job0']); data=fread(pathBruker,precision); fclose(pathBruker); %read fid
data=data(1:2:end)+1i*data(2:2:end); %complex data
data=reshape(data,Np,sum(AverageList),NR);

if NR~=1; data=squeeze(sum(data,3));end

kspace=zeros(Np,AcqMatrix(1)*AcqMatrix(2)); k=0;
for p=1:length(AverageList)
    for l=1:AverageList(p)
        k=k+1;
        kspace(:,p)=kspace(:,p)+data(:,k);
    end
end
kspace=reshape(kspace,Np,AcqMatrix(2),AcqMatrix(1));
kspace=conj(kspace);

rawdata.Nx=AcqMatrix(2);
rawdata.Ny=AcqMatrix(1);
rawdata.Np=Np;
rawdata.SW=SW;
rawdata.kspace=kspace;

