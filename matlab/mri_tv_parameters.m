function [methods,audit] = mri_tv_parameters(kind,candidates,gridName)
%MRI_TV_PARAMETERS Predeclared grids independent of image/test results.
% publication48: 48 configurations per method. Includes all balanced12 settings.
% MRI/RI: 4 sigma x 3 alpha x 2 primal-step fractions x 2 relaxation fractions.
% CV/CP: 4 sigma x 12 primal-step fractions; no inertia and rho=1.
if nargin<2, candidates=true; end
if nargin<3, gridName='publication48'; end
kind=upper(char(kind)); gridName=lower(char(gridName));
assert(any(strcmp(gridName,{'publication48','balanced12','paper12'})), ...
    'MRI:Grid','Unknown grid.');
methods=struct('kind',{},'alpha',{},'rho',{},'tau',{},'sigma',{},'id',{});
if ~candidates
    switch kind
        case {'MRI','RI'}, methods=make(kind,.2,.3,.2,.3,1);
        case 'CV', methods=make(kind,0,1,.9/(.5+8*.3),.3,1);
        case 'CP', methods=make(kind,0,1,.9/(9*.3),.3,1);
        otherwise, error('MRI:Method','Unknown method.');
    end
else
    id=0;
    switch kind
        case {'MRI','RI'}
            if strcmp(kind,'MRI'), alphas=[0 .3 .9]; else, alphas=[0 .1 .3]; end
            if strcmp(gridName,'publication48')
                for a=alphas
                    for sig=[.03 .1 .3 1]
                        root=2/(8*sig+sqrt((8*sig)^2+4));
                        for stepFraction=[.6 .9]
                            tau=stepFraction*root;
                            for fraction=[.9 .99]
                                rho=fraction*rho_bound(kind,a,tau,sig);
                                id=id+1; methods(id)=make(kind,a,rho,tau,sig,id);
                            end
                        end
                    end
                end
            elseif strcmp(gridName,'balanced12')
                for a=alphas
                    for sig=[.03 .1 .3 1]
                        root=2/(8*sig+sqrt((8*sig)^2+4));
                        tau=.9*root;rho=.9*rho_bound(kind,a,tau,sig);
                        id=id+1;methods(id)=make(kind,a,rho,tau,sig,id);
                    end
                end
            else
                sig=.3;root=2/(8*sig+sqrt((8*sig)^2+4));
                for a=alphas
                    for fraction=[.5 .9]
                        for stepFraction=[.6 .9]
                            tau=stepFraction*root;rho=fraction*rho_bound(kind,a,tau,sig);
                            id=id+1;methods(id)=make(kind,a,rho,tau,sig,id);
                        end
                    end
                end
            end
        case {'CV','CP'}
            if strcmp(gridName,'publication48')
                fractions=[.3 .4 .5 .6 .7 .8 .85 .9 .93 .95 .97 .99];
            else
                fractions=[.3 .6 .9];
            end
            for sig=[.03 .1 .3 1]
                for fraction=fractions
                    if strcmp(kind,'CV'),tau=fraction/(.5+8*sig);
                    else,tau=fraction/(9*sig);end
                    id=id+1;methods(id)=make(kind,0,1,tau,sig,id);
                end
            end
        otherwise,error('MRI:Method','Unknown method.');
    end
end
audit=cell(1,numel(methods));
for j=1:numel(methods),audit{j}=mri_tv_check_parameters(methods(j));end
end
function r=rho_bound(kind,a,tau,sig)
if strcmp(kind,'MRI'),r=1/(1+a);return;end
ell=tau/sqrt(1-8*tau*sig);epsilon=(1-ell^2)/100;
psi=(2-epsilon+2*ell)/(1+ell^2+2*ell);
r=psi*(1-a)^2/(2*a^2-a+1);
end
function m=make(kind,a,r,t,s,id)
m=struct('kind',kind,'alpha',a,'rho',r,'tau',t,'sigma',s,'id',id);
end
