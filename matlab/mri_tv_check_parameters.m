function a = mri_tv_check_parameters(m)
%MRI_TV_CHECK_PARAMETERS Conservative sufficient conditions for fixed blur/TV.
flds={'kind','alpha','rho','tau','sigma'};
assert(isstruct(m) && isscalar(m) && all(isfield(m,flds)), ...
    'MRI:InvalidParameters','Invalid parameter structure.');
validateattributes(m.tau,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(m.sigma,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(m.alpha,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(m.rho,{'numeric'},{'real','finite','scalar','positive'});
a=struct();
switch m.kind
    case {'MRI','RI'}
        a.operatorBudget=8*m.tau*m.sigma+m.tau^2;
        assert(a.operatorBudget<1,'MRI:Inadmissible','PD kernel bound failed.');
        a.ell=m.tau/sqrt(1-8*m.tau*m.sigma);
        if strcmp(m.kind,'MRI')
            a.rhoMax=1/(1+m.alpha);
            a.parameterBudget=m.rho*(1+m.alpha);
            assert(a.parameterBudget<1,'MRI:Inadmissible','MRI inertia/relaxation bound failed.');
        else
            assert(m.alpha<1,'MRI:Inadmissible','RI requires alpha<1.');
            a.epsilon=(1-a.ell^2)/100;
            a.psi=(2-a.epsilon+2*a.ell)/(1+a.ell^2+2*a.ell);
            a.rhoMax=a.psi*(1-m.alpha)^2/(2*m.alpha^2-m.alpha+1);
            assert(m.rho<a.rhoMax,'MRI:Inadmissible','RI inertia/relaxation bound failed.');
        end
    case 'CV'
        a.operatorBudget=m.tau*(.5+8*m.sigma);
        assert(a.operatorBudget<1 && m.rho==1 && m.alpha==0, ...
            'MRI:Inadmissible','Condat-Vu bound/parameters failed.');
    case 'CP'
        a.operatorBudget=9*m.tau*m.sigma;
        assert(a.operatorBudget<1 && m.rho==1 && m.alpha==0, ...
            'MRI:Inadmissible','Chambolle-Pock bound/parameters failed.');
    otherwise, error('MRI:Method','Unknown method: %s',m.kind);
end
end
