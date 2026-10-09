function d = mri_tv_diagnostics(x,u,b,lambdaTV,op,cfg)
%MRI_TV_DIAGNOSTICS Common KKT residual and a genuine primal-dual gap.
% The KKT residual is tested at the RAW iterate, not at a clipped replacement.
% Objective/PSNR are reported at xf=projection_[0,255](x), and feasibility
% violations of the raw iterate are separately recorded and tested.
% Gap uses vf=K*xf-b and uf=projection_[-lambdaTV,lambdaTV](u) for ALL solvers.
box=@(z) min(255,max(0,z));
dual=@(z) min(lambdaTV,max(-lambdaTV,z));
nrm=@(z) norm(z(:));
rx=op.Kt(op.K(x)-b)+op.Gt(u);
rx=(x-box(x-cfg.diagnosticTau*rx))/cfg.diagnosticTau;
ru=(u-dual(u+cfg.diagnosticSigma*op.G(x)))/cfg.diagnosticSigma;
d.kkt=max(nrm(rx)/max(1,nrm(x)),nrm(ru)/max(1,nrm(u)));
xf=box(x); uf=dual(u);
v=op.K(xf)-b; gx=op.G(xf);
d.objective=0.5*sum(v(:).^2)+lambdaTV*sum(abs(gx(:)));
s=-op.Kt(v)-op.Gt(uf);
% Stable evaluation of F(xf)-dual(v,uf) as two nonnegative Fenchel gaps.
tvgap=sum(lambdaTV*abs(gx(:))-gx(:).*uf(:));
boxgap=sum((255-xf(:)).*max(s(:),0)-xf(:).*min(s(:),0));
gap=tvgap+boxgap;
if gap < -1e-10*max(1,d.objective)
    error('Negative primal-dual gap: operator/projection implementation error.');
end
d.gap=max(0,gap);
d.dualLower=d.objective-d.gap;
d.relativeGap=d.gap/max(1,abs(d.objective));
d.primalViolation=max(abs(x(:)-xf(:)))/255;
d.dualViolation=max(abs(u(:)-uf(:)))/max(1,lambdaTV);
d.feas=max(d.primalViolation,d.dualViolation);
d.success=isfinite(d.kkt) && isfinite(d.relativeGap) && ...
    d.kkt<=cfg.tolKKT && d.relativeGap<=cfg.tolGap && d.feas<=cfg.tolFeas;
end
