function [s,aux] = mri_tv_step(s,b,lambdaTV,op,m,collectResidual)
%MRI_TV_STEP One update. MRI/RI use the PROVED FBF primal-dual realization.
% CV is Condat-Vu, not incorrectly relabelled as MRI-NFB.
if nargin<6, collectResidual=false; end
box=@(z) min(255,max(0,z));
dual=@(z) min(lambdaTV,max(-lambdaTV,z));
aux.inclusionResidual=NaN;
xold=s.x; uold=s.u;
switch m.kind
    case {'MRI','RI'}
        p=s.x+m.alpha*(s.x-s.xprev);
        q=s.u+m.alpha*(s.u-s.uprev);
        fp=op.Kt(op.K(p)-b);
        xt=box(p-m.tau*(fp+op.Gt(q)));
        fxt=op.Kt(op.K(xt)-b);
        w=xt-m.tau*(fxt-fp);
        ut=dual(q+m.sigma*op.G(xt+w-p));
        if strcmp(m.kind,'MRI')
            s.x=(1-m.rho)*xold+m.rho*w;
            s.u=(1-m.rho)*uold+m.rho*ut;
        else
            s.x=(1-m.rho)*p+m.rho*w;
            s.u=(1-m.rho)*q+m.rho*ut;
        end
        if collectResidual
            % r=M(p,q)-M(xt,ut)=S*((p,q)-(w,ut))/tau at (xt,ut).
            ex=p-w; eu=q-ut; gex=op.G(ex);
            metric2=sum(ex(:).^2)-2*m.tau*sum(gex(:).*eu(:))+ ...
                (m.tau/m.sigma)*sum(eu(:).^2);
            assert(metric2>=-1e-9*max(1,sum(ex(:).^2)), ...
                'Negative product-metric norm.');
            aux.inclusionResidual=sqrt(max(0,metric2))/m.tau;
            aux.xback=xt; aux.uback=ut;
        end
    case 'CV'
        xt=box(s.x-m.tau*(op.Kt(op.K(s.x)-b)+op.Gt(s.u)));
        ut=dual(s.u+m.sigma*op.G(2*xt-s.x));
        s.x=xt; s.u=ut;
    case 'CP'
        s.v=(s.v+m.sigma*(op.K(s.xbar)-b))/(1+m.sigma);
        s.u=dual(s.u+m.sigma*op.G(s.xbar));
        s.x=box(s.x-m.tau*(op.Kt(s.v)+op.Gt(s.u)));
        s.xbar=2*s.x-xold;
    otherwise
        error('Unknown solver kind: %s',m.kind);
end
s.xprev=xold; s.uprev=uold;
end
