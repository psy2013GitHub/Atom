
function [nodes,matrix,sig]=NBS4strucNet(hc_data,mdd_data,NPerm,Zthresh,Pthresh)
hc_rawData=importdata(hc_data);
mdd_rawData=importdata(mdd_data);
%## init ##%
Roinames=hc_rawData.colheaders; nroi=length(Roinames);
hc_n=size(hc_rawData.data);   hc_n=hc_n(1);
mdd_n=size(mdd_rawData.data); mdd_n=mdd_n(1);
rawData=cat(1,hc_rawData.data,mdd_rawData.data);
%## permutation loop##%
for n =0:NPerm
    fprintf('.')
    if ~n
        norder=1:(hc_n+mdd_n);
    else
        norder=randperm(hc_n+mdd_n);
    end
    hc_data=rawData(norder(1:hc_n),:);
    mdd_data=rawData(norder(hc_n+1:end),:);
    hc_M=pearsonM(hc_data);
    mdd_M=pearsonM(mdd_data);
    if ~n
        actual_hc_M=hc_M;
        actual_mdd_M=mdd_M;
    end
    outM = r_test(hc_n,mdd_n,hc_M,mdd_M,1);
    outM(outM<Zthresh)=0; outM(outM>=Zthresh)=1;
    [comps,comp_sizes] = get_components(outM);
    if diff(size(comp_sizes))>0, comp_sizes=comp_sizes'; end
    if ~n
        actual_comps=comps;
        actual_comp_sizes=comp_sizes;
        glCnt=zeros(1,length(actual_comp_sizes));
    else
        comp_sizes_compa=bsxfun(@le,actual_comp_sizes',comp_sizes);
        glCnt = glCnt+logical(sum(comp_sizes_compa,1));
    end
end
glCnt=glCnt/NPerm;
sig_comps=find(glCnt<=Pthresh);
fprintf(['\n',num2str(length(sig_comps)),' Components Found\n']);
N_sig_comps=length(sig_comps);
nodes=cell(N_sig_comps,1); matrix=cell(N_sig_comps,1); sig=cell(N_sig_comps,1);
for ii=1:length(sig_comps)
    nodes{ii}={Roinames{find(actual_comps==sig_comps(ii))}};
    matrix{ii}=cat(acutal_hc_M(nods,nods),acutal_mdd_M(nods,nods));
    sig{ii}=glCnt(sig_comps(ii));
end
disp();

end

function rM = pearsonM(rawdata)
[nr,nc]=size(rawdata);
rawdata=rawdata-ones(nr,1)*mean(rawdata);
rawdata=rawdata./(ones(nr,1)*std(rawdata));
rM=rawdata'*rawdata/(nr-1);
rM=rM-diag(diag(rM));
end

function ZM = fisherZ(rM)
ZM=0.5 * (log(1+rM) - log(1-rM));
end

function [z,p] = r_test(n1,n2,r1,r2,twotail_option)
z1=fisherZ(r1);
z2=fisherZ(r2);
se_diff_r=(1/(n1-3)+1/(n2-3))^0.5;
z=abs((z1-z2)./se_diff_r);
p=1-(erf(z/(2^0.5))+1)/2;
if twotail_option
    p=p*2;
end
end

function [comps,comp_sizes] = get_components(adj)
if size(adj,1)~=size(adj,2)
    error('this adjacency matrix is not square');
end

if ~any(adj-triu(adj))
  adj = adj | adj';
end

%if main diagonal of adj do not contain all ones, i.e. autoloops
if sum(diag(adj))~=size(adj,1)
    adj = adj|speye(size(adj));
end

%Dulmage-Mendelsohn decomposition
[useless1,p,useless2,r] = dmperm(adj);

comp_sizes = diff(r);
num_comps = numel(comp_sizes);
comps = zeros(1,size(adj,1)); 
comps(r(1:num_comps)) = ones(1,num_comps); 
comps = cumsum(comps); 
comps(p) = comps; 
end