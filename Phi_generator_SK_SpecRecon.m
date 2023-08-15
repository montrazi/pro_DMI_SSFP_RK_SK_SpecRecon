function Phi = Phi_generator_SK_SpecRecon(timeexp,K)

beta_ini    = 0.05;
beta_end    = 0.5; %mu

KTrans      = 1;

kep_glu_ini = 0.015;
kep_glu_end = 0.1;

kep_pyr_ini = 0.015;
kep_pyr_end = 0.1;

kep_lac_ini = 0.015;
kep_lac_end = 0.1;

kep_wat_ini = 0.000005;
kep_wat_end = 0.0005;

if isrow(timeexp); timeexp=timeexp.'; end
t = timeexp;

%Glucose
plato = logspace(log10(0.2),log10(1),40);
backgr = linspace(0,0.1,5);
betalist = logspace(log10(beta_ini),log10(beta_end),5);
kepGlulist = logspace(log10(kep_glu_ini),log10(kep_glu_end),5);
Cglu = zeros(length(backgr),length(betalist),length(plato),length(kepGlulist),length(t));
for iback = 1:length(backgr)
    for ibeta = 1:length(betalist)
        for ip = 1:length(plato)
            for iGlu = 1:length(kepGlulist)
                Cp = exp(-betalist(ibeta)*t);
                kepGlu = kepGlulist(iGlu);
                Cglu(iback,ibeta,ip,iGlu,:) = TrapzKetyFun(Cp,[KTrans kepGlu],t)';
                Cglu(iback,ibeta,ip,iGlu,:) = Cglu(iback,ibeta,ip,iGlu,:)/squeeze(max(Cglu(iback,ibeta,ip,iGlu,:),[],'all'));
                Cglu(iback,ibeta,ip,iGlu,:) = min(Cglu(iback,ibeta,ip,iGlu,:),plato(ip));
                Cglu(iback,ibeta,ip,iGlu,:) = Cglu(iback,ibeta,ip,iGlu,:)/squeeze(max(Cglu(iback,ibeta,ip,iGlu,:),[],'all'));
                Cglu(iback,ibeta,ip,iGlu,:) = ( Cglu(iback,ibeta,ip,iGlu,:) + backgr(iback)/(1-backgr(iback)) ) / ( 1 + backgr(iback)/(1-backgr(iback)) );
            end
        end
    end
end
CgluFinal = reshape(Cglu,[size(Cglu,1)*size(Cglu,2)*size(Cglu,3)*size(Cglu,4),size(Cglu,5)]);
CgluFinal(:,1) = 0;

%Pyruvate
Cglu = reshape(Cglu(1,:,:,:,:),[1*size(Cglu,2)*size(Cglu,3)*size(Cglu,4),size(Cglu,5)]);
Cglulist = 1:round(size(Cglu,1)/100):size(Cglu,1);
kepPyrlist = logspace(log10(kep_pyr_ini),log10(kep_pyr_end),5);
Cpyr = zeros(size(Cglulist,1),length(kepPyrlist),length(t));
for iGlu = 1:length(Cglulist)
    for  iPyr = 1:length(kepPyrlist)
        Cg = Cglu(Cglulist(iGlu),:).';
        kepPyr = kepPyrlist(iPyr);
        Cpyr(iGlu,iPyr,:) = TrapzKetyFun(Cg,[KTrans kepPyr],t)';
        Cpyr(iGlu,iPyr,:) = Cpyr(iGlu,iPyr,:)/squeeze(max(Cpyr(iGlu,iPyr,:),[],'all'));
    end
end
Cpyr = reshape(Cpyr,[size(Cpyr,1)*size(Cpyr,2),size(Cpyr,3)]);

%Lactate
backgr = linspace(0,0.3,10);
Cpyrlist = 1:round(size(Cpyr,1)/250):size(Cpyr,1);
kepLaclist = logspace(log10(kep_lac_ini),log10(kep_lac_end),5);
Clac = zeros(length(backgr),size(Cpyrlist,1),length(kepLaclist),length(t));
for iback = 1:length(backgr)
    for iPyr = 1:length(Cpyrlist)
        for  iLac = 1:length(kepLaclist)
            Cp = Cpyr(Cpyrlist(iPyr),:).';
            kepLac = kepLaclist(iLac);
            Clac(iback,iPyr,iLac,:) = TrapzKetyFun(Cp,[KTrans kepLac],t)';
            Clac(iback,iPyr,iLac,:) = Clac(iback,iPyr,iLac,:)/squeeze(max(Clac(iback,iPyr,iLac,:),[],'all'));
            Clac(iback,iPyr,iLac,:) = ( Clac(iback,iPyr,iLac,:) + backgr(iback)/(1-backgr(iback)) ) / ( 1 + backgr(iback)/(1-backgr(iback)) );
        end
    end
end
Clac = reshape(Clac,[size(Clac,1)*size(Clac,2)*size(Clac,3),size(Clac,4)]);

%Water
backgr = linspace(0,0.5,10);
Cpyrlist = 1:round(size(Cpyr,1)/250):size(Cpyr,1);
kepWatlist = logspace(log10(kep_wat_ini),log10(kep_wat_end),5);
Cwat = zeros(length(backgr),size(Cpyrlist,1),length(kepWatlist),length(t));
for iback = 1:length(backgr)
    for iPyr = 1:length(Cpyrlist)
        for  iWat = 1:length(kepWatlist)
            Cp = Cpyr(Cpyrlist(iPyr),:).';
            kepWat = kepWatlist(iWat);
            Cwat(iback,iPyr,iWat,:) = TrapzKetyFun(Cp,[KTrans kepWat],t)';
            Cwat(iback,iPyr,iWat,:) = Cwat(iback,iPyr,iWat,:)/squeeze(max(Cwat(iback,iPyr,iWat,:),[],'all'));
            Cwat(iback,iPyr,iWat,:) = ( Cwat(iback,iPyr,iWat,:) + backgr(iback)/(1-backgr(iback)) ) / ( 1 + backgr(iback)/(1-backgr(iback)) );
        end
    end
end
Cwat = reshape(Cwat,[size(Cwat,1)*size(Cwat,2)*size(Cwat,3),size(Cwat,4)]);

Cglu=CgluFinal;

% figure; plot(t,Cwat,'Color','black');
figure; plot(t,Cwat(1:end,:));
xlabel('time (min)'); ylabel('C_{HDO}'); title('HDO concentration in the tissues'); set(gca,'FontSize',14);
% set(gca,'visible','off');
X = Cwat';
[U,~,~] = svd(X,'econ');
Phi_wat = U(:,1:K);
% colorphi={'#0072BD' '#D95319' '#EDB120' '#7E2F8E'};
% figure; for p=1:4; subplot(4,1,p); plot(t,Phi_wat(:,p),'Color',char(colorphi(p))); set(gca,'FontSize',14); end
figure; plot(t,Phi_wat); xlabel('time (min)'); ylabel('\phi_{HDO}'); set(gca,'FontSize',14);
Z = Phi_wat * Phi_wat' * X;
err = norm(X(:) - Z(:)) / norm(X(:));
% fprintf('Relative norm of error: %.6f\n', err);

% figure; plot(t,Cglu,'Color','black');
figure; plot(t,Cglu(:,:));
xlabel('time (min)'); ylabel('C_{Glu}'); title('Glucose concentration in the tissues'); set(gca,'FontSize',14);
% set(gca,'visible','off');
X = Cglu';
[U,~,~] = svd(X,'econ');
Phi_glu = U(:,1:K);
figure; plot(t,Phi_glu); xlabel('time (min)'); ylabel('\phi_{Glu}'); set(gca,'FontSize',14);
Z = Phi_glu * Phi_glu' * X;
err = norm(X(:) - Z(:)) / norm(X(:));
% fprintf('Relative norm of error: %.6f\n', err);

% figure; plot(t,Clac,'Color','black');
figure; plot(t,Clac(:,:));
xlabel('time (min)'); ylabel('C_{Lac}'); title('Lactate concentration in the tissues'); set(gca,'FontSize',14);
% set(gca,'visible','off');
X = Clac';
[U,~,~] = svd(X,'econ');
Phi_lac = U(:,1:K);
% colorphi={'#0072BD' '#D95319' '#EDB120' '#7E2F8E'};
% figure; for p=1:4; subplot(4,1,p); plot(t,Phi_lac(:,p),'Color',char(colorphi(p))); set(gca,'FontSize',14); end
figure; plot(t,Phi_lac); xlabel('time (min)'); ylabel('\phi_{Lac}'); set(gca,'FontSize',14);
Z = Phi_lac * Phi_lac' * X;
err = norm(X(:) - Z(:)) / norm(X(:));
% fprintf('Relative norm of error: %.6f\n', err);

Phi(:,:,1) = Phi_wat;
Phi(:,:,2) = Phi_glu;
Phi(:,:,3) = Phi_lac;

% print('YourEPSFile','-dpdf','-vector');
end
