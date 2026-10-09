function op = mri_tv_operators(sz, backend)
%MRI_TV_OPERATORS Exact Neumann gradient/adjoint and reflective 3x3 blur.
% K = B_h X B_w with B symmetric, row-stochastic, three-point averaging.
% The FFT backend uses an even 2h-by-2w extension. It is NOT a preconditioner
% and does NOT replace Neumann conditions by periodic boundary conditions.
if nargin < 2, backend = 'spatial'; end
validateattributes(sz, {'numeric'}, {'vector','numel',2,'integer','>=',2});
h = sz(1); w = sz(2); backend = char(backend);
if ~any(strcmp(backend, {'spatial','fft_reflect'})), error('Invalid backend.'); end
nK=0; nKt=0; nG=0; nGt=0; nFFT=0;
if strcmp(backend,'fft_reflect')
    vr=(1+2*cos(2*pi*(0:2*h-1)'/(2*h)))/3;
    vc=(1+2*cos(2*pi*(0:2*w-1)/(2*w)))/3;
    transfer=vr*vc;
else
    transfer=[];
end
op.K=@applyK; op.Kt=@applyKt; op.G=@applyG; op.Gt=@applyGt;
op.counts=@getCounts; op.reset=@resetCounts;
op.normK2=1; op.normG2=8; op.size=[h w]; op.backend=backend;
    function y=blur(x)
        if ~isequal(size(x),[h w]), error('K: image size mismatch.'); end
        if strcmp(backend,'spatial')
            y=(x(:,[1 1:w-1])+x+x(:,[2:w w]))/3;
            y=(y([1 1:h-1],:)+y+y([2:h h],:))/3;
        else
            e=[x,fliplr(x);flipud(x),flipud(fliplr(x))];
            e=real(ifft2(fft2(e).*transfer));
            y=e(1:h,1:w); nFFT=nFFT+2;
        end
    end
    function y=applyK(x), nK=nK+1; y=blur(x); end
    function y=applyKt(x), nKt=nKt+1; y=blur(x); end
    function u=applyG(x)
        nG=nG+1;
        if ~isequal(size(x),[h w]), error('G: image size mismatch.'); end
        u=zeros(h,w,2);
        u(:,1:w-1,1)=x(:,2:w)-x(:,1:w-1);
        u(1:h-1,:,2)=x(2:h,:)-x(1:h-1,:);
    end
    function x=applyGt(u)
        nGt=nGt+1;
        if ~isequal(size(u),[h w 2]), error('Gt: field size mismatch.'); end
        % Ignore the dummy final column/row: G is exactly zero there.
        x=zeros(h,w);
        x(:,1:w-1)=x(:,1:w-1)-u(:,1:w-1,1);
        x(:,2:w)=x(:,2:w)+u(:,1:w-1,1);
        x(1:h-1,:)=x(1:h-1,:)-u(1:h-1,:,2);
        x(2:h,:)=x(2:h,:)+u(1:h-1,:,2);
    end
    function c=getCounts()
        c=struct('K',nK,'Kt',nKt,'G',nG,'Gt',nGt,'FFT',nFFT);
    end
    function resetCounts(), nK=0; nKt=0; nG=0; nGt=0; nFFT=0; end
end
