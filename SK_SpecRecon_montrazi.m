function outParams=SK_SpecRecon_montrazi(imDataParams,algoParams,forceB0,Phi_matrix)

Nimages=imDataParams.images; %imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
TEs=imDataParams.TE; if isrow(TEs); TEs=TEs.'; end %imDataParams.TE: echo times (in seconds)
Ny=size(Nimages,1); Nx=size(Nimages,2); Nz=size(Nimages,3);

GyromagneticRatio=42.576; % MHz/T
LarmorFreq=imDataParams.FieldStrength*GyromagneticRatio; %imDataParams.FieldStrength: (in Tesla)

Nspec=numel(algoParams.species);
wfreq=zeros(Nspec,1);

mul=-imDataParams.PrecessionIsClockwise;
for n=1:Nspec
     wfreq(n,1)=2*pi*LarmorFreq*mul*algoParams.species(n).frequency*algoParams.species(n).relAmps(:); %ppm to frequency (in rad/s)
end

%%%% Make the matrix A %%%%
NTEs=size(TEs,1); A=zeros(NTEs,Nspec);
for j=1:NTEs
    A(j,:)=exp(1i*wfreq*TEs(j,1))';
end
A=kron(eye(Nz),A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Subspace for kinetic dimension %%%%
K=size(Phi_matrix,2);
Phi=zeros(Nspec*Nz,Nspec*K);
j=1;
for jj=1:Nz
    for n=1:Nspec
        Phi(j,(n-1)*K+1:(n-1)*K+K)=squeeze(Phi_matrix(jj,:,n));
        j=j+1;
    end
end
% assignin('base','Phi',Phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Make inversion matrix %%%%
Mat=(A*Phi)'*(A*Phi);
MaxCond=1000;
[U,S,V]=svd(Mat,0);
S=diag(S);
idx1=find(S>S/MaxCond);
idx2=find(S<=S/MaxCond);
CondNum=S(1)/min(S(idx1));
S(idx1)=1./S(idx1);
S(idx2)=0;
InvMat=V*diag(S)*U';
if isempty(S(idx2))
    fprintf('Condition number for inversion = %f\n',CondNum);
else
    fprintf('Ill-conditioned inversion matrix (condition number limited to %f)\n',CondNum);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nimages=permute(Nimages,[1 2 5 3 4]);
Nimages=reshape(Nimages,Ny,Nx,Nz*NTEs);
SpiMat=spiral_index(Nimages);
pspec=zeros(Ny,Nx,Nz,Nspec);
TEs=repmat(TEs,Nz,1);
B0map=zeros(Ny,Nx,Nz);
R2map=zeros(Ny,Nx,Nz);
if strcmp(forceB0,'off')
for n=1:Ny*Nx
    if rem(n,100)==0; fprintf('Loop %d\n out of %d\n',n,Ny*Nx); end
    s0=squeeze(Nimages(SpiMat(n,1),SpiMat(n,2),:));
    if s0~=0
    phi=zeros(Nz,1);
    if SpiMat(n,1)>=3 && SpiMat(n,2)>=3 && SpiMat(n,1)<=Ny-2 && SpiMat(n,2)<=Nx-2
        phi=B0map(SpiMat(n,1)-2,SpiMat(n,2)-2,:)+B0map(SpiMat(n,1)-2,SpiMat(n,2)-1,:)+B0map(SpiMat(n,1)-2,SpiMat(n,2)+0,:)+B0map(SpiMat(n,1)-2,SpiMat(n,2)+1,:)+B0map(SpiMat(n,1)-2,SpiMat(n,2)+2,:)+...
            B0map(SpiMat(n,1)-1,SpiMat(n,2)-2,:)+B0map(SpiMat(n,1)-1,SpiMat(n,2)-1,:)+B0map(SpiMat(n,1)-1,SpiMat(n,2)+0,:)+B0map(SpiMat(n,1)-1,SpiMat(n,2)+1,:)+B0map(SpiMat(n,1)-1,SpiMat(n,2)+2,:)+...
            B0map(SpiMat(n,1)+0,SpiMat(n,2)-2,:)+B0map(SpiMat(n,1)+0,SpiMat(n,2)-1,:)+B0map(SpiMat(n,1)+0,SpiMat(n,2)+0,:)+B0map(SpiMat(n,1)+0,SpiMat(n,2)+1,:)+B0map(SpiMat(n,1)+0,SpiMat(n,2)+2,:)+...
            B0map(SpiMat(n,1)+1,SpiMat(n,2)-2,:)+B0map(SpiMat(n,1)+1,SpiMat(n,2)-1,:)+B0map(SpiMat(n,1)+1,SpiMat(n,2)+0,:)+B0map(SpiMat(n,1)+1,SpiMat(n,2)+1,:)+B0map(SpiMat(n,1)+1,SpiMat(n,2)+2,:)+...
            B0map(SpiMat(n,1)+2,SpiMat(n,2)-2,:)+B0map(SpiMat(n,1)+2,SpiMat(n,2)-1,:)+B0map(SpiMat(n,1)+2,SpiMat(n,2)+0,:)+B0map(SpiMat(n,1)+2,SpiMat(n,2)+1,:)+B0map(SpiMat(n,1)+2,SpiMat(n,2)+2,:);
        phi=squeeze(phi)/25;
    end
    windowSize=3;
    b=(1/windowSize)*ones(1,windowSize);
    a=1;
    phi=filter(b,a,phi);
    for rep=1:10
        s=s0.*exp(-1i*2*pi*kron(phi,ones(NTEs,1)).*TEs);
%         sigma=((A*Phi)'*(A*Phi))\(A*Phi)'*s;
        sigma=InvMat*(A*Phi)'*s;
        rho=Phi*sigma;
        er=s-A*rho;
            
        iT=2*pi*1i*diag(TEs);
        rhotilde=zeros(Nspec*Nz,Nz);
        for k=1:Nz
            rhotilde((k-1)*Nspec+1:k*Nspec,k)=rho((k-1)*Nspec+1:k*Nspec,1);
        end
        B=iT*A*rhotilde;

        %%%% Make inversion matrix %%%%
        MatB=B'*B;
        [U,S,V]=svd(MatB,0);
        S=diag(S);
        idx1=find(S>S/MaxCond);
        idx2=find(S<=S/MaxCond);
        S(idx1)=1./S(idx1);
        S(idx2)=0;
        InvMatB=V*diag(S)*U';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         Dphi=(B'*B)\B'*er;
        Dphi=InvMatB*B'*er;
        phi=phi+real(Dphi);
    end
    end
    B0map(SpiMat(n,1),SpiMat(n,2),:)=phi;
end
end

B0map=imboxfilt3(B0map);

errormap=zeros(size(B0map));
for n=1:Ny
    for m=1:Nx
        phi=squeeze(B0map(n,m,:));
        s=squeeze(Nimages(n,m,:));
        s=s.*exp(-1i*2*pi*kron(phi,ones(NTEs,1)).*TEs);
%         sigma=((A*Phi)'*(A*Phi))\(A*Phi)'*s;
        sigma=InvMat*(A*Phi)'*s;
        rho=Phi*sigma;
        pspec(n,m,:,:)=reshape(rho,Nspec,Nz).';
        errormap(n,m,:)=sum(reshape(s-A*rho,NTEs,Nz),1);
    end
end

for n=1:Nspec
    outParams.species(n).amps(:,:,:,1)=pspec(:,:,:,n);
end

outParams.B0map=B0map;
outParams.r2starmap=R2map;
outParams.errormap=errormap;

end

