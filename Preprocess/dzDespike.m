
function data = dzDespike(data,TR,Options)

c1 = 2.5;
c2 = 3;
if exist('Options','var')
    if isfield(Options,'c1'), c1=Options.c1; end
    if isfield(Options,'c2'), c2=Options.c2; end
end

nTp=size(data,1); nCol=size(data,2);

[lestimates]            = dzRegress(data,[ones(nTp,1) (-1:2/(nTp-1):1)']);
[qestimates,  modelq]   = icatb_myquadfun(data,TR);
[splestimates,  models] = icatb_mysplinefun(data,TR);


ylfit =  lestimates(1,:) + lestimates(2,:)*(-1:2/(nTp-1):1)';
yqfit = icatb_getQuadFit(qestimates,length(tc),TR);
ysfit = icatb_getSplineFit(splestimates,length(tc),TR);

err = [icatb_gfit2(tc,ylfit,'1') icatb_gfit2(tc,yqfit,'1') icatb_gfit2(tc,ysfit,'1')];

[mnerr mnID] = min(err);

if mnID == 1
    yfit =  ylfit;
elseif mnID == 2
    yfit = yqfit;
else
    yfit = ysfit;
end

res = tc - yfit;
mad_res = median(abs(res - median(res))); % median absolute deviation of residuals
sigma = mad_res* sqrt(pi/2);
s = res/sigma;
tc_out = s;

ind = find(abs(s) > c1);
tc_out(ind) = sign(s(ind))*(c1+((c2-c1)*tanh((abs(s(ind))-c1)/(c2-c1))));

tc_out = yfit + tc_out*sigma;

end

function [estimates,  model] = icatb_myquadfun(V,TR)
%TR = 2;
nt = length(V);

%numP = floor(nt/30);
t = 0:TR:(nt-1)*TR;
t = t(:);
[Bhat] = icatb_regress(V,[ones(nt,1) t]);
start_point = [0 Bhat(2) Bhat(1)];
model = @myspfun;
options = optimset('MaxFunEvals', 100000, 'Display', 'off');
estimates = fminsearch(model, start_point, options);

    function [err yfun] = myspfun(params)
        % despike AFNI method; estimate a smoothish curve

        x0 = params;

        % smoothish quadratic fit
        yfun = x0(1)*t.^2 + x0(2)*t + x0(3);
        %         for ii = 1:numP
        %             yfun = yfun + x0(2+ii) * sin(2*pi*ii*t/(nt*TR)) + x0(2+ii) *cos(2*pi*ii*t/(nt*TR));
        %         end

        err = sum((V - yfun).^2);

    end

end

function [estimates,  model] = icatb_mysplinefun(V,TR)
%TR = 2;
nt = length(V);

numP = floor(nt/30);
t = 0:TR:(nt-1)*TR;
t = t(:);
start_point = rand(numP+2, 1);
model = @myspfun;
options = optimset('MaxFunEvals', 100000, 'Display', 'off');
estimates = fminsearch(model, start_point, options);

    function [err yfun] = myspfun(params)
        % despike AFNI method; estimate a smoothish curve

        x0 = params;

        % smoothish spline fit
        yfun = x0(1)*t + x0(2)*t.^2;
        for ii = 1:numP
            yfun = yfun + x0(2+ii) * sin(2*pi*ii*t/(nt*TR)) + x0(2+ii) *cos(2*pi*ii*t/(nt*TR));
        end

        err = sum((V - yfun).^2);
    end
end

