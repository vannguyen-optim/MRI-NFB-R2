function checks = test_mri_tv(runConvergence)
%TEST_MRI_TV Unit checks; NOT experimental evidence for bird.tif.
% test_mri_tv(false): algebra/adjoints/FFT/parameter tests only.
% test_mri_tv(true): also solve a small synthetic problem using all 4 methods.
if nargin<1, runConvergence=true; end
oldrng=rng; cleanup=onCleanup(@()rng(oldrng)); %#ok<NASGU>
rng(701,'twister');
a=mri_tv_operators([7 9],'spatial'); f=mri_tv_operators([7 9],'fft_reflect');
x=randn(7,9); y=randn(7,9); u=randn(7,9,2);
ip=@(s,t) sum(s(:).*t(:));
rel=@(s,t) abs(s-t)/max([1 abs(s) abs(t)]);
checks.adjointG=rel(ip(a.G(x),u),ip(x,a.Gt(u)));
checks.adjointK=rel(ip(a.K(x),y),ip(x,a.Kt(y)));
checks.adjointKFFT=rel(ip(f.K(x),y),ip(x,f.Kt(y)));
checks.fftAgreement=norm(a.K(x)-f.K(x),'fro')/max(1,norm(x,'fro'));
checks.constants=norm(a.K(ones(7,9))-ones(7,9),'fro');
assert(max([checks.adjointG checks.adjointK checks.adjointKFFT checks.fftAgreement checks.constants])<1e-12);
% Exact KKT solution: a constant observation and zero dual variable.
cfg=mri_tv_config('smoke'); b=100*ones(7,9); lam=3;
d=mri_tv_diagnostics(b,zeros(7,9,2),b,lam,a,cfg);
checks.exactKKT=d.kkt; checks.exactGap=d.gap;
assert(d.kkt<1e-12 && d.gap<1e-10 && d.objective<1e-10);
for grid={'publication48','balanced12','paper12'}
    for kind={'MRI','RI','CV','CP'}
        ms=mri_tv_parameters(kind{1},true,grid{1});
        expected=12;if strcmp(grid{1},'publication48'),expected=48;end
        assert(numel(ms)==expected);
        for k=1:numel(ms), mri_tv_check_parameters(ms(k)); end
    end
end
% General Fenchel gap equals the independently evaluated dual formula.
xr=-10+280*rand(7,9); ur=8*randn(7,9,2); br=100+10*randn(7,9);
dr=mri_tv_diagnostics(xr,ur,br,lam,a,cfg);
xf=min(255,max(0,xr)); uf=min(lam,max(-lam,ur)); vv=a.K(xf)-br;
ss=-a.Kt(vv)-a.Gt(uf);
dualValue=-.5*sum(vv(:).^2)-sum(br(:).*vv(:))-255*sum(max(ss(:),0));
checks.dualFormula=abs(dr.gap-(dr.objective-dualValue))/max(1,abs(dr.gap));
assert(checks.dualFormula<1e-12);
% A cap need not be divisible by checkEvery; diagnostics must match the return.
short=cfg; short.maxIter=7; short.checkEvery=10; short.tolKKT=1e-14; short.tolGap=1e-14;
mm=mri_tv_parameters('MRI',false);
[xx,uu,out]=mri_tv_solve(br,lam,a,mm,short);
last=mri_tv_diagnostics(xx,uu,br,lam,a,short);
assert(out.history(end,1)==out.iter && isequaln(out.final,last));
checks.finalDiagnostics=true;
% Dense product-space check with NONLINEAR monotone Lipschitz D.
% This is independent of the imaging gradient used by mri_tv_step.
tau=.15; sig=.2; L=[.3 -.2 .1; .1 .4 -.3]; n=3; m=2;
D=@(z) .5*z+.1*sin(z);
S=[eye(n),-tau*L';-tau*L,(tau/sig)*eye(m)];
M=@(z) [z(1:n)/tau-D(z(1:n))-L'*z(n+1:end); ...
    -L*z(1:n)+tau*L*D(z(1:n))+z(n+1:end)/sig];
z=[.4;-.8;1.6;.2;-.4]; p=z(1:n); q=z(n+1:end);
xt=min(1,max(0,p-tau*(D(p)+L'*q)));
w=xt-tau*(D(xt)-D(p)); ut=min(lam,max(-lam,q+sig*L*(xt+w-p)));
zt=[xt;ut]; zhat=z-tau*(S\(M(z)-M(zt)));
checks.productCorrection=norm(zhat-[w;ut]); assert(checks.productCorrection<1e-12);
r=M(z)-M(zt); normalx=r(1:n)-D(xt)-L'*ut;
normalu=r(n+1:end)+L*xt;
checks.primalNormal=norm(xt-min(1,max(0,xt+normalx)));
checks.dualNormal=norm(ut-min(lam,max(-lam,ut+normalu)));
assert(checks.primalNormal<1e-12 && checks.dualNormal<1e-12);
checks.metricResidual=abs(sqrt(r'*(S\r))-sqrt((z-zhat)'*S*(z-zhat))/tau);
assert(checks.metricResidual<1e-12);
% Alpha=0: MRI and RI are the same algorithm at matched parameters.
b=100+10*randn(7,9); s.x=b; s.xprev=b+randn(7,9);
s.u=randn(7,9,2); s.uprev=randn(7,9,2); s.v=zeros(7,9); s.xbar=b;
mi=mri_tv_parameters('MRI',false); mi.alpha=0;
mj=mi; mj.kind='RI';
si=mri_tv_step(s,b,lam,a,mi); sj=mri_tv_step(s,b,lam,a,mj);
checks.zeroInertia=norm(si.x-sj.x,'fro')+norm(si.u(:)-sj.u(:));
assert(checks.zeroInertia<1e-12);
% Corrected MRI/RI must agree between two implementations of the SAME K.
si=mri_tv_step(s,b,lam,a,mi); sf=mri_tv_step(s,b,lam,f,mi);
checks.fftStep=norm(si.x-sf.x,'fro')+norm(si.u(:)-sf.u(:));
assert(checks.fftStep<1e-10);
% Boundary illustration: T=-I, A=0, B=I, mu=2, M=I/2.
x=1; xp=-.5; alpha=2; rho=.2;
for k=1:300
    y=x+alpha*(x-xp); xn=(1-rho)*x-rho*y; xp=x; x=xn;
end
checks.boundaryLinear=abs(x); assert(abs(x)<1e-10);
if runConvergence
    [rr,cc]=ndgrid(1:12,1:14);
    truth=80+60*(rr>4 & rr<10 & cc>3 & cc<11);
    a=mri_tv_operators(size(truth),'spatial'); b=a.K(truth)+4*randn(size(truth));
    cfg=mri_tv_config('benchmark'); cfg.maxIter=30000;
    values=zeros(1,4); kinds={'MRI','RI','CV','CP'};
    for j=1:4
        mm=mri_tv_parameters(kinds{j},false);
        [~,~,res]=mri_tv_solve(b,lam,a,mm,cfg);
        assert(res.success,'Synthetic convergence test did not meet the tolerances.');
        values(j)=res.final.objective;
        checks.([kinds{j} '_iterations'])=res.iter;
    end
    checks.objectiveAgreement=(max(values)-min(values))/max(1,min(values));
    assert(checks.objectiveAgreement<5e-4,'Solvers disagree beyond gap tolerance.');
end
checks.passed=true;
fprintf('All requested mri_tv checks passed. These tests are not publication results.\n');
end
