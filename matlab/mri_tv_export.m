function mri_tv_export(records,selected,calibration,manifest,cfg,finalize)
%MRI_TV_EXPORT Raw output only. Publication table/figure export is separate.
% The finalize argument is retained for API compatibility.
write_raw(fullfile(cfg.outDir,'raw_results.csv'),records);
write_raw(fullfile(cfg.outDir,'calibration_results.csv'),calibration);
fid=openout(fullfile(cfg.outDir,'selected_parameters.csv')); guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Image,Method,Candidate,Alpha,Rho,Tau,Sigma\n');
for i=1:size(selected,1)
    for j=1:size(selected,2)
        m=selected{i,j}; if isempty(m), continue; end
        fprintf(fid,'%d,%s,%d,%.17g,%.17g,%.17g,%.17g\n',i,m.kind,m.id,m.alpha,m.rho,m.tau,m.sigma);
    end
end
if finalize && isempty(manifest), error('MRI:Manifest','Missing manifest.'); end
end
function fid=openout(path)
fid=fopen(path,'w'); if fid<0, error('MRI:WriteFailed','Cannot write %s',path); end
end
function write_raw(path,records)
fid=openout(path); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Panel,Image,Seed,Backend,Method,Alpha,Rho,Tau,Sigma,LambdaTV,Iterations,ElapsedSeconds,Success,Status,ObjectivePerPixel,PSNR,KKT,RelativeGap,Feasibility,PrimalViolation,DualViolation,DualLowerBound,K,Kt,G,Gt,FFT\n');
for k=1:numel(records)
    z=records{k}; m=z.method; d=z.final; c=z.counts;
    fprintf(fid,'%s,%d,%d,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d\n', ...
        z.panel,z.image,z.seed,z.backend,m.kind,m.alpha,m.rho,m.tau,m.sigma, ...
        z.lambdaTV,z.iter,z.elapsed,z.success,z.status,d.objective/z.pixels,z.psnr, ...
        d.kkt,d.relativeGap,d.feas,d.primalViolation,d.dualViolation,d.dualLower, ...
        c.K,c.Kt,c.G,c.Gt,c.FFT);
end
end
