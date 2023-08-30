function data=OpenBrukerImage(pathfile)

visu_parsBruker=fopen([pathfile '\visu_pars']); p=1; infoBruker=[]; infoBruker{1}=fgetl(visu_parsBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(visu_parsBruker); end
fclose(visu_parsBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1)
    if contains(infoBruker{p},'##$VisuCoreWordType='); VisuCoreWordType=erase(infoBruker{p},'##$VisuCoreWordType='); end
    if contains(infoBruker{p},'##$VisuCoreByteOrder='); VisuCoreByteOrder=erase(infoBruker{p},'##$VisuCoreByteOrder='); end
    if contains(infoBruker{p},'##$VisuCoreDataSlope=')
        VisuCoreDataSlope=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; VisuCoreDataSlope=[VisuCoreDataSlope str2num(infoBruker{pp})]; pp=pp+1; end
    end
    if contains(infoBruker{p},'##$VisuCoreSize=')
        VisuCoreSize=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; VisuCoreSize=[VisuCoreSize str2num(infoBruker{pp})]; pp=pp+1; end
    end
    if contains(infoBruker{p},'##$VisuCoreDataSlope=')
        VisuCoreDataOffs=[]; pp=p+1;
        while contains(infoBruker{pp},'##$')~=1; VisuCoreDataOffs=[VisuCoreDataOffs str2num(infoBruker{pp})]; pp=pp+1; end
    end
end
numDataHighDim=prod(VisuCoreSize(2:end));
switch VisuCoreWordType
    case '_32BIT_SGN_INT'; precision='int32';
    case '_16BIT_SGN_INT'; precision='int16';
    case '_8BIT_UNSGN_INT'; precision='uint8';
    case '_32BIT_FLOAT'; precision='single';
end
switch VisuCoreByteOrder
    case 'littleEndian'; endian='l';
    case 'bigEndian'; endian='b';
end

recoBruker=fopen([pathfile '\reco']); p=1; infoBruker=[]; infoBruker{1}=fgetl(recoBruker);
while ischar(infoBruker{p}) p=p+1; infoBruker{p}=fgetl(recoBruker); end
fclose(recoBruker); infoBruker=infoBruker(1:p-1)';
for p=1:size(infoBruker,1)
    if contains(infoBruker{p},'##$RECO_ft_size'); sizmat=[str2num(infoBruker{p+1}) str2num(infoBruker{p+2})]; end
end

pathBruker=fopen([pathfile '\2dseq']); data=fread(pathBruker,precision,0,endian); fclose(pathBruker); %read 2dseq
data=reshape(data,sizmat(1),sizmat(2),length(data)/sizmat(1)/sizmat(2));

slope=repmat(VisuCoreDataSlope',[VisuCoreSize(1)*numDataHighDim,1]);
slope=reshape(slope,size(data));
data=data.*slope;
offset=repmat(VisuCoreDataOffs',[VisuCoreSize(1)*numDataHighDim,1]);
offset=reshape(offset,size(data));
data=data+offset;