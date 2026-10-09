function [x,u,out] = mri_tv_solve(b,lambdaTV,op,m,cfg)
%MRI_TV_SOLVE Time includes ALL stopping diagnostics and their operator calls.
% Data generation and shared setup are excluded; output I/O is excluded.
mri_tv_check_parameters(m);
s.x=b; s.xprev=b; s.u=zeros([size(b) 2]); s.uprev=s.u;
s.v=zeros(size(b)); s.xbar=b;
bNorm=max(1,norm(b(:)));
op.reset(); tstart=tic;
d=mri_tv_diagnostics(s.x,s.u,b,lambdaTV,op,cfg);
hist=zeros(ceil(cfg.maxIter/cfg.checkEvery)+2,7); nh=1;
hist(nh,:)=[0,toc(tstart),d.objective,d.kkt,d.relativeGap,d.feas,NaN];
success=d.success; iter=0; status='maxiter';
if success, status='converged'; end
for k=1:cfg.maxIter
    if success, break; end
    collect=(mod(k,cfg.checkEvery)==0 || k==cfg.maxIter);
    [s,aux]=mri_tv_step(s,b,lambdaTV,op,m,collect);
    iter=k;
    if any(~isfinite(s.x(:))) || any(~isfinite(s.u(:))) || ...
            any(~isfinite(s.v(:))) || any(~isfinite(s.xbar(:))) || ...
            norm(s.x(:))>1e12*bNorm
        status='numerical_failure'; break;
    end
    if collect
        d=mri_tv_diagnostics(s.x,s.u,b,lambdaTV,op,cfg);
        nh=nh+1;
        hist(nh,:)=[k,toc(tstart),d.objective,d.kkt,d.relativeGap, ...
            d.feas,aux.inclusionResidual];
        success=d.success;
        if success, status='converged'; end
    end
end
out.elapsed=toc(tstart); out.iter=iter; out.success=success;
out.status=status; out.history=hist(1:nh,:); out.final=d;
out.counts=op.counts(); out.method=m;
out.historyColumns={'iteration','elapsedSeconds','objective','kkt','relativeGap','feasibility','inclusionResidual'};
if strcmp(status,'numerical_failure')
    out.final.kkt=Inf; out.final.relativeGap=Inf;
    out.final.feas=Inf; out.final.success=false;
    out.final.objective=NaN; out.final.gap=Inf; out.final.dualLower=NaN;
    out.final.primalViolation=Inf; out.final.dualViolation=Inf;
    out.success=false;
    if out.history(end,1)~=iter
        out.history(end+1,:)=[iter,out.elapsed,NaN,Inf,Inf,Inf,NaN];
    end
end
x=s.x; u=s.u;
end
