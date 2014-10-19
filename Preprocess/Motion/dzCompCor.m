
function [WM,CSF]=dzCompCor(Functionals,Masks)

nSubjs=length(Functionals);
% cCompCor/-> wm/csf mask erosion
for subj=1:nSubjs,
    ERODE=1; % erosing level for white/csf masks (voxels) (set to 0 for no erosion)
    if ERODE>0,
        Vmask{nroi}=conn_prepend('e',CONN_x.Setup.rois.files{nsub}{nroi}{1});
        for m=1:length(Masks{subj})
            [nill,nill,ext]=fileparts(Masks{subj}{m});
            switch(ext),
                case {'.img','.nii','.hdr'},
                    V0=spm_vol(Masks{subj}{m}); % mask
                    X0=spm_read_vols(V0);
                    idx1=find(X0(:)>.5);
                    [idxx,idxy,idxz]=ind2sub(size(X0),idx1);
                    idxt=find(idxx>ERODE&idxx<size(X0,1)+1-ERODE&idxy>ERODE&idxy<size(X0,2)+1-ERODE&idxz>ERODE&idxz<size(X0,3)+1-ERODE);
                    for n1=1:length(idxt),
                        if (sum(sum(sum(X0(idxx(idxt(n1))+(-ERODE:ERODE),idxy(idxt(n1))+(-ERODE:ERODE),idxz(idxt(n1))+(-ERODE:ERODE))<.5,3),2),1))>1, 
                            idxt(n1)=0; 
                        end; 
                    end
                    idxt=idxt(idxt>0);
                    idx1=idx1(idxt);
                    X1=zeros(size(X0));X1(idx1)=1;
                    V0.fname=conn_prepend('e',Masks{subj}{m}{1});spm_write_vol(V0,X1);
                otherwise,
                    Vmask{nroi}=CONN_x.Setup.rois.files{nsub}{nroi}{1}; % fix me
            end
        end
    end
end

% rex 里面提取主成分
if strcmpi(params.summary_measure,'eigenvariate'), % first-step: compute covariance structure
    data=zeros(length(params.VF),length(params.VF));
    datamean=zeros(1,size(XYZmm,2));
    dataM=zeros(length(params.VF),1);
    dataN=0;
    if ~silence, hft=waitbar(0,['Precomputing covariance structure']); set(hft,'color','w'); end
    for n1=1:1e3:size(XYZmm,2),
        idx=n1:min(size(XYZmm,2),n1-1+1e3);
        temp1=zeros(length(params.VF),length(idx));
        for i=1:length(params.VF),
            XYZ = iM{i}*[XYZmm(:,idx); ones(1,length(idx))];
            temp1(i,:) = spm_sample_vol(params.VF(i),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
        end
        idxnan=find(isnan(temp1));
        if ~isempty(idxnan),
            datamean(idx)=sum(~isnan(temp1),1);
            dataN=dataN+sum(datamean(idx)>0);
            temp1(idxnan)=0;
            datamean(idx)=sum(temp1,1)./max(eps,datamean(idx));
            [idxnani,idxnanj]=ind2sub(size(temp1),idxnan);
            temp1(idxnan)=datamean(idx(idxnanj));
        else datamean(idx)=mean(temp1,1); dataN=dataN+length(idx); end
        if ~silence, waitbar(n1/size(XYZmm,2),hft); end
        temp2=temp1-repmat(mean(temp1,1),[size(temp1,1),1]);
        data=data+temp2*temp2';
        dataM=dataM+sum(temp2,2);
    end
    %             data=data/dataN;
    %             dataM=dataM/dataN;
    %             data=data-dataM*dataM';
    if ~silence, close(hft); end
    cov0=ones(size(data,1),1);
    if isfield(params,'pca')&&~params.pca, dataM=[]; end
    if isfield(params,'covariates')&&~isempty(params.covariates),
        cov1=[dataM,detrend(params.covariates,'constant')];
    else
        cov1=dataM;
    end
    proj=eye(size(data,1))-[cov1,0*cov0]*pinv([cov1,cov0]); % removes covariates keeping scale unchanged
    data=proj*data*proj';
    [q1,q2,nill]=svd(data);
    if isfield(params,'pca')&&~params.pca, temp=min(size(q1,2),params.dims(min(r,length(params.dims))));
    else temp=min(size(q1,2),params.dims(min(r,length(params.dims)))-1);
    end
    data=q1(:,1:temp)*diag(sqrt(diag(q2(1:temp,1:temp))));
    basis=zeros(size(XYZmm,2),temp);
end