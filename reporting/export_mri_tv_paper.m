function report = export_mri_tv_paper(resultsDir,makePlots,figuresOnly)
%EXPORT_MRI_TV_PAPER Read completed run data and export numerical tables/figures.
% No retuning, no new solver execution, and no pre-written performance claim.
% File names match the manually fillable Section 4. Captions remain in TeX.
% Graphics patch: same white-background style as the synthetic experiment.
% Each line chart is saved as vector PDF and 600-DPI PNG.
% Native restored-image PNG files are kept; PDF companions are added.
% Set figuresOnly=true to leave all saved tables and numerical data untouched.
% Examples:
%   export_mri_tv_paper(resultsDir,true)       % tables and figures
%   export_mri_tv_paper(resultsDir,true,true)  % figures only
if nargin<2, makePlots=true; end
if nargin<3, figuresOnly=false; end
file=fullfile(resultsDir,'benchmark.mat');
assert(exist(file,'file')==2,'MRI:MissingResults','Missing benchmark.mat: %s',resultsDir);
z=load(file,'R'); R=z.R; cfg=R.config; mf=R.manifest;
assert(R.complete,'MRI:Incomplete','Experiment is incomplete. Resume the run first.');
records=R.records; nImages=size(R.selected,1); kinds=cfg.methodKinds;
for nm={'tables','figures'}
    p=fullfile(resultsDir,nm{1}); if ~exist(p,'dir'), mkdir(p); end
end
if figuresOnly
    assert(makePlots,'MRI:PlotOption','figuresOnly requires makePlots=true.');
    for im=1:nImages
        plot_histories(R,im,resultsDir);
        plot_sensitivity(R,im,resultsDir);
        export_image_pdfs(im,resultsDir);
    end
    report=struct('outputDirectory',resultsDir,'images',nImages, ...
        'solverRerun',false,'figuresOnly',true);
    fprintf('TV figures exported as PDF and PNG: %s\n',fullfile(resultsDir,'figures'));
    return;
end
issues={};
for im=1:nImages
    for j=1:numel(kinds)
        rr=take(records,im,'benchmark',kinds{j},''); seeds=values(rr,'seed');
        assert(numel(unique(seeds))==numel(seeds),'MRI:Duplicates','Duplicate test records.');
        if ~isequal(sort(seeds(:)),sort((cfg.seedBase+(1:mf.nRuns))'))
            issues{end+1}=sprintf('Image %d %s: test records incomplete.',im,kinds{j}); %#ok<AGROW>
        end
    end
end
if any(cellfun(@(r)strcmp(r.status,'numerical_failure'),records))
    issues{end+1}='Numerical failures occurred. Inspect raw_results.csv; do not omit them.';
end
if strcmp(cfg.profile,'publication')
    if mf.nRuns~=20 || ~isequal(cfg.calibrationSeeds,[1001 1002 1003]) || cfg.seedBase~=2026
        issues{end+1}='The saved seed protocol differs from the reference publication protocol.';
    end
    if cfg.maxIter~=20000 || cfg.tolKKT~=1e-4 || cfg.tolGap~=1e-4 || cfg.tolFeas~=1e-8
        issues{end+1}='Caps/tolerances differ from the reference manuscript.';
    end
end
fid=openout(fullfile(resultsDir,'protocol_audit.txt')); guard=onCleanup(@()fclose(fid));
fprintf(fid,'Profile: %s\nGrid: %s\nImages: %d\nTest realizations/image: %d\n',cfg.profile,cfg.gridName,nImages,mf.nRuns);
fprintf(fid,'MATLAB: %s\nPlatform: %s\n',mf.version,mf.computer);
if isempty(issues), fprintf(fid,'No completeness or numerical-failure issue found by this audit.\n'); end
for k=1:numel(issues), fprintf(fid,'%s\n',issues{k}); end
clear guard;
benchRows={}; accuracyRows={}; sensitivityRows={}; backendRows={};
for im=1:nImages
    tag=sprintf('i%02d',im);
    fid=openout(fullfile(resultsDir,'tables',['paper_parameters_' tag '.tex'])); guard=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n','\begin{tabular}{lccccc}\hline');
    row(fid,{'Method','$\tau$','$\sigma$','$\alpha$','$\rho$','Calibration success'});
    fprintf(fid,'%s\n','\hline');
    for j=1:numel(kinds)
        m=choice(R.selected,im,kinds{j}); cal=take(R.calibration,im,'calibration',kinds{j},'');
        cal=cal(cellfun(@(r)same_method(r.method,m),cal));
        a=number(m.alpha,5); r=number(m.rho,6); cs='--';
        if cfg.tune, cs=count_success(cal); end
        if any(strcmp(kinds{j},{'CV','CP'})), a='--'; r='--'; end
        row(fid,{method_name(kinds{j}),number(m.tau,8),number(m.sigma,8),a,r,cs});
    end
    fprintf(fid,'%s\n','\hline\end{tabular}'); clear guard;
    fid=openout(fullfile(resultsDir,'tables',['paper_benchmark_' tag '.tex'])); guard=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n','\begin{tabular}{lccccc}\hline');
    row(fid,{'Method','Success','Iterations','Time (s)','Operator calls','PSNR (dB)'});
    fprintf(fid,'%s\n','\hline');
    for j=1:numel(kinds)
        rr=take(records,im,'benchmark',kinds{j},''); if isempty(rr), continue; end
        nk=cellfun(@(r)r.counts.K+r.counts.Kt,rr); ng=cellfun(@(r)r.counts.G+r.counts.Gt,rr);
        work=sprintf('%.1f / %.1f',mean(nk),mean(ng));
        row(fid,{method_name(kinds{j}),count_success(rr),pm(values(rr,'iter'),1),pm(values(rr,'elapsed'),3),work,pm(values(rr,'psnr'),4)});
        benchRows(end+1,:)={im,kinds{j},numel(rr),sum(values(rr,'success')), ...
            sum(cellfun(@(r)strcmp(r.status,'maxiter'),rr)), ...
            sum(cellfun(@(r)strcmp(r.status,'numerical_failure'),rr)), ...
            mean(values(rr,'iter')),std(values(rr,'iter'),0),mean(values(rr,'elapsed')), ...
            std(values(rr,'elapsed'),0),mean(nk),mean(ng), ...
            mean(values(rr,'psnr')),std(values(rr,'psnr'),0)}; %#ok<AGROW>
    end
    fprintf(fid,'%s\n','\hline\end{tabular}'); clear guard;
    fid=openout(fullfile(resultsDir,'tables',['paper_accuracy_' tag '.tex'])); guard=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n','\begin{tabular}{lcccc}\hline');
    row(fid,{'Method','Objective/pixel','KKT residual','Relative gap','Feasibility'});
    fprintf(fid,'%s\n','\hline');
    for j=1:numel(kinds)
        rr=take(records,im,'benchmark',kinds{j},''); if isempty(rr), continue; end
        obj=cellfun(@(r)r.final.objective/r.pixels,rr); kkt=median_all(finalvals(rr,'kkt'));
        gap=median_all(finalvals(rr,'relativeGap')); feas=median_all(finalvals(rr,'feas'));
        row(fid,{method_name(kinds{j}),pm(obj,4),sci(kkt),sci(gap),sci(feas)});
        accuracyRows(end+1,:)={im,kinds{j},mean(obj),std(obj,0),kkt,gap,feas}; %#ok<AGROW>
    end
    fprintf(fid,'%s\n','\hline\end{tabular}'); clear guard;
    fid=openout(fullfile(resultsDir,'tables',['paper_sensitivity_' tag '.tex'])); guard=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n','\begin{tabular}{ccccccc}\hline');
    row(fid,{'$\rho$','$\alpha$','MRI / RI success','Pairs','$\Delta$PSNR','Iteration ratio','Time ratio'});
    fprintf(fid,'%s\n','\hline');
    for rho=cfg.sensitivityRho
        for a=cfg.sensitivityAlpha
            aa=setting(records,im,'MRI',a,rho); bb=setting(records,im,'RI',a,rho);
            if isempty(aa) || isempty(bb), continue; end
            delta=[]; ir=[]; tr=[];
            for t=1:numel(aa)
                idx=find(cellfun(@(r)r.seed==aa{t}.seed,bb),1);
                if isempty(idx), continue; end
                u=aa{t}; v=bb{idx};
                if ~u.success || ~v.success || u.iter<=0 || u.elapsed<=0, continue; end
                delta(end+1)=u.psnr-v.psnr; ir(end+1)=v.iter/u.iter; tr(end+1)=v.elapsed/u.elapsed; %#ok<AGROW>
            end
            row(fid,{number(rho,3),number(a,3),[count_success(aa) ' / ' count_success(bb)], ...
                sprintf('%d',numel(tr)),pm(delta,4),pm(ir,3),pm(tr,3)});
            sensitivityRows(end+1,:)={im,rho,a,sum(values(aa,'success')),sum(values(bb,'success')), ...
                numel(aa),numel(bb),numel(tr),avg(delta),sdval(delta),avg(ir),sdval(ir),avg(tr),sdval(tr)}; %#ok<AGROW>
        end
    end
    fprintf(fid,'%s\n','\hline\end{tabular}'); clear guard;
    fid=openout(fullfile(resultsDir,'tables',['paper_backend_' tag '.tex'])); guard=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n','\begin{tabular}{lcccc}\hline');
    row(fid,{'Implementation','Success','Iterations','Time (s)','Relative gap'});
    fprintf(fid,'%s\n','\hline');
    for be={'spatial','fft_reflect'}
        rr=take(records,im,'backend','MRI',be{1}); if isempty(rr), continue; end
        if strcmp(be{1},'spatial'), label='Direct reflection'; else, label='FFT with even reflection'; end
        row(fid,{label,count_success(rr),pm(values(rr,'iter'),1),pm(values(rr,'elapsed'),3),sci(median_all(finalvals(rr,'relativeGap')))});
        backendRows(end+1,:)={im,be{1},sum(values(rr,'success')),numel(rr),mean(values(rr,'iter')), ...
            std(values(rr,'iter'),0),mean(values(rr,'elapsed')),std(values(rr,'elapsed'),0),median_all(finalvals(rr,'relativeGap'))}; %#ok<AGROW>
    end
    fprintf(fid,'%s\n','\hline\end{tabular}'); clear guard;
    if makePlots
        plot_histories(R,im,resultsDir);
        plot_sensitivity(R,im,resultsDir);
        export_image_pdfs(im,resultsDir);
    end
end
write_cells(fullfile(resultsDir,'summary_benchmark.csv'), ...
    {'Image','Method','Runs','Success','Capped','NumericalFailures','IterationsMean','IterationsSD', ...
    'TimeMeanSeconds','TimeSDSeconds','K_KtMeanCalls','G_GtMeanCalls','PSNRMean','PSNRSD'},benchRows);
write_cells(fullfile(resultsDir,'summary_accuracy.csv'), ...
    {'Image','Method','ObjectivePerPixelMean','ObjectivePerPixelSD','KKTMedian','RelativeGapMedian','FeasibilityMedian'},accuracyRows);
write_cells(fullfile(resultsDir,'summary_sensitivity.csv'), ...
    {'Image','Rho','Alpha','MRISuccess','RISuccess','MRIRuns','RIRuns','Pairs','DeltaPSNRMean','DeltaPSNRSD', ...
    'IterationRatioMean','IterationRatioSD','TimeRatioMean','TimeRatioSD'},sensitivityRows);
write_cells(fullfile(resultsDir,'summary_backend.csv'), ...
    {'Image','Implementation','Success','Runs','IterationsMean','IterationsSD','TimeMeanSeconds','TimeSDSeconds','RelativeGapMedian'},backendRows);
report=struct('outputDirectory',resultsDir,'images',nImages,'protocolIssues',{issues},'solverRerun',false);
fprintf('CSV summaries and LaTeX fragments written to %s\n',resultsDir);
end
function y=avg(v), if isempty(v), y=NaN; else, y=mean(v); end, end
function y=sdval(v), if isempty(v), y=NaN; else, y=std(v,0); end, end
function write_cells(path,header,rows)
fid=openout(path); guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',strjoin(header,','));
for i=1:size(rows,1)
    out=cell(1,size(rows,2));
    for j=1:size(rows,2)
        v=rows{i,j};
        if isnumeric(v)||islogical(v), out{j}=sprintf('%.17g',v);
        else, out{j}=['"' strrep(char(v),'"','""') '"']; end
    end
    fprintf(fid,'%s\n',strjoin(out,','));
end
end
function plot_histories(R,im,root)
% Plot the saved first-test histories, not newly computed solver output.
kinds={'MRI','RI','CV','CP'};
seed=R.config.seedBase+1;
colors=curve_colors();
lineStyles={'-','--','-.',':'};
markerShapes={'o','s','^','d'};
for col=[4 5]
    [f,ax]=paper_axes('log');
    guard=onCleanup(@()close_figure(f));
    lineHandles=gobjects(0,1);
    leg={};
    lastTime=0;
    for j=1:numel(kinds)
        rr=take(R.records,im,'benchmark',kinds{j},'');
        ix=find(cellfun(@(z)z.seed==seed,rr),1);
        if isempty(ix), continue; end
        rec=rr{ix}; h=rec.history;
        if isempty(h), continue; end
        % Same time and diagnostic columns as in the original exporter.
        valid=isfinite(h(:,2)) & isfinite(h(:,col)) & h(:,col)>=0;
        if ~any(valid), continue; end
        yy=max(h(valid,col),realmin);
        xx=h(valid,2);
        lastTime=max(lastTime,max(xx));
        % All samples are plotted; only the marker positions are staggered.
        firstMarker=min(3+2*(j-1),numel(xx));
        markerEvery=max(12,ceil(numel(xx)/8));
        markerIndices=firstMarker:markerEvery:numel(xx);
        lineHandles(end+1,1)=plot(ax,xx,yy, ...
            'Color',colors(j,:), ...
            'LineStyle',lineStyles{j}, ...
            'Marker',markerShapes{j}, ...
            'MarkerIndices',markerIndices, ...
            'MarkerSize',4, ...
            'MarkerEdgeColor',colors(j,:), ...
            'MarkerFaceColor',colors(j,:), ...
            'LineWidth',1.3); %#ok<AGROW>
        leg{end+1}=sprintf('%s (%s)',plain_name(kinds{j}),rec.status); %#ok<AGROW>
    end
    if isempty(lineHandles)
        clear guard;
        continue;
    end
    if col==4
        yLabel='Common KKT residual';
        tag='kkt_time'; tol=R.config.tolKKT;
    else
        yLabel='Relative primal-dual gap';
        tag='gap_time'; tol=R.config.tolGap;
    end
    label_axes(ax,'Elapsed time (s)',yLabel);
    % Use the actual saved times. No hard-coded limits from Experiment 1.
    if lastTime>0, xlim(ax,[0 lastTime]); else, xlim(ax,[0 1]); end
    lim=xlim(ax);
    plot(ax,lim,[tol tol],'Color','k','LineStyle','--', ...
        'LineWidth',1.0,'HandleVisibility','off');
    paper_legend(ax,lineHandles,leg,'southwest','none');
    save_figure(f,fullfile(root,'figures',sprintf('paper_i%02d_%s.pdf',im,tag)));
    clear guard;
end
end

function plot_sensitivity(R,im,root)
% Retain the original aggregation and the original missing-setting gaps.
if isempty(take(R.records,im,'sensitivity','','')), return; end
kinds={'MRI','RI'};
colors=curve_colors();
markers={'o','s','^','d','v','>'};
methodStyles={'-','--'};
nRho=numel(R.config.sensitivityRho);
for isTime=[false true]
    [f,ax]=paper_axes('linear');
    guard=onCleanup(@()close_figure(f));
    leg={}; lineHandles=gobjects(0,1);
    for ir=1:nRho
        rho=R.config.sensitivityRho(ir);
        thisColor=colors(1+mod(ir-1,size(colors,1)),:);
        for j=1:numel(kinds)
            aa=R.config.sensitivityAlpha;
            yy=NaN(size(aa));
            for k=1:numel(aa)
                rr=setting(R.records,im,kinds{j},aa(k),rho);
                if isempty(rr), continue; end
                if isTime
                    yy(k)=median_all(values(rr,'elapsed'));
                else
                    yy(k)=sum(cellfun(@(z)z.success,rr))/numel(rr);
                end
            end
            finiteIndices=find(isfinite(yy));
            if isempty(finiteIndices), continue; end
            styleIndex=2*(ir-1)+j;
            markerIndex=1+mod(styleIndex-1,numel(markers));
            % One or more actual data markers per curve; no x/y offsets.
            start=min(styleIndex,numel(finiteIndices));
            markerIndices=finiteIndices(start:max(1,2*nRho):end);
            if j==1, markerFace=thisColor; else, markerFace='w'; end
            lineHandles(end+1,1)=plot(ax,aa,yy, ...
                'Color',thisColor,'LineStyle',methodStyles{j}, ...
                'Marker',markers{markerIndex}, ...
                'MarkerIndices',markerIndices, ...
                'MarkerSize',4,'MarkerEdgeColor',thisColor, ...
                'MarkerFaceColor',markerFace,'LineWidth',1.3); %#ok<AGROW>
            leg{end+1}=sprintf('%s, \\rho = %.3g',plain_name(kinds{j}),rho); %#ok<AGROW>
        end
    end
    if isTime
        yLabel='Median elapsed time (s)';
        tag='sensitivity_time'; location='best';
    else
        yLabel='Success proportion';
        tag='sensitivity_success'; location='southwest';
        ylim(ax,[-0.03 1.03]);
        yticks(ax,0:0.2:1);
    end
    label_axes(ax,'Inertial parameter \alpha',yLabel);
    % Use the configured domain, including inadmissible settings left as NaN.
    alphas=R.config.sensitivityAlpha;
    if numel(alphas)>1 && min(alphas)<max(alphas)
        xlim(ax,[min(alphas) max(alphas)]);
    end
    if ~isempty(lineHandles)
        paper_legend(ax,lineHandles,leg,location,'tex');
    end
    save_figure(f,fullfile(root,'figures',sprintf('paper_i%02d_%s.pdf',im,tag)));
    clear guard;
end
end

function [f,ax]=paper_axes(yScale)
% Match run_mri_synthetic.m: white canvas, black text, Arial, light grid.
f=figure('Visible','off','Units','inches', ...
    'Position',[1 1 6.2 4.2],'Color','w', ...
    'Name','TV restoration','NumberTitle','off');
if isprop(f,'Theme'), theme(f,'light'); end
set(f,'Color','w');
ax=axes('Parent',f,'Units','normalized', ...
    'Position',[0.14 0.15 0.83 0.80], ...
    'Color','w','XColor','k','YColor','k', ...
    'XScale','linear','YScale',yScale, ...
    'FontName','Arial','FontSize',10,'FontWeight','normal', ...
    'LineWidth',0.75,'TickDir','out', ...
    'XMinorTick','off','YMinorTick','off');
% Set YScale BEFORE hold on (important for logarithmic histories).
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'GridColor','k','GridAlpha',0.2,'GridLineStyle','-', ...
    'XMinorGrid','off','YMinorGrid','off','Layer','top');
end

function label_axes(ax,xLabel,yLabel)
xlabel(ax,xLabel,'FontName','Arial','FontSize',11, ...
    'FontWeight','normal','Color','k','Interpreter','tex');
ylabel(ax,yLabel,'FontName','Arial','FontSize',11, ...
    'FontWeight','normal','Color','k','Interpreter','tex');
end

function paper_legend(ax,lineHandles,labels,location,interpreter)
legend(ax,lineHandles,labels,'Location',location, ...
    'Interpreter',interpreter,'FontName','Arial','FontSize',10, ...
    'Color','w','TextColor','k','EdgeColor','k', ...
    'Box','on','AutoUpdate','off');
end

function colors=curve_colors()
% Same four base colours as the final synthetic figure.
colors=[31 119 180; 255 127 14; 44 160 44; 214 39 40]/255;
end

function save_figure(f,file)
% Both files come from the SAME figure. No modification of numerical data.
[folder,name]=fileparts(file);
pdfFile=fullfile(folder,[name '.pdf']);
pngFile=fullfile(folder,[name '.png']);
drawnow;
exportgraphics(f,pdfFile, ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(f,pngFile, ...
    'Resolution',600,'BackgroundColor','white');
fprintf('Saved figures:\n%s\n%s\n',pdfFile,pngFile);
end

function export_image_pdfs(im,root)
% Keep the native-resolution grayscale PNGs generated by the original code.
% Add PDF companions using these exact saved PNG pixels. No contrast change,
% inversion, normalization, new observation, or solver call is performed.
% Photographs remain raster images embedded in their PDF (not vector curves).
figDir=fullfile(root,'figures');
labels={'truth','observed','MRI','RI','CV','CP'};
for j=1:numel(labels)
    name=sprintf('i%02d_%s',im,labels{j});
    pngFile=fullfile(figDir,[name '.png']);
    if exist(pngFile,'file')~=2, continue; end
    [pixels,map]=imread(pngFile);
    h=size(pixels,1); w=size(pixels,2);
    if h==0 || w==0, continue; end
    figureSize=4.2*[w h]/max(w,h);
    f=figure('Visible','off','Units','inches', ...
        'Position',[1 1 figureSize],'Color','w','NumberTitle','off');
    guard=onCleanup(@()close_figure(f));
    if isprop(f,'Theme'), theme(f,'light'); end
    set(f,'Color','w');
    ax=axes('Parent',f,'Units','normalized','Position',[0 0 1 1], ...
        'Color','w','XColor','k','YColor','k');
    if ~isempty(map)
        image(ax,pixels,'CDataMapping','direct');
        colormap(ax,map);
    elseif ndims(pixels)==2
        % The publication pipeline writes uint8 PNGs: intensity 0..255.
        assert(isa(pixels,'uint8'),'MRI:PlotImageType', ...
            'Expected the native uint8 grayscale output: %s',pngFile);
        image(ax,pixels,'CDataMapping','direct');
        colormap(ax,gray(256));
    else
        image(ax,pixels);
    end
    axis(ax,'image'); axis(ax,'off');
    set(ax,'YDir','reverse','XLim',[0.5 w+0.5],'YLim',[0.5 h+0.5]);
    drawnow;
    pdfFile=fullfile(figDir,[name '.pdf']);
    exportgraphics(ax,pdfFile,'ContentType','vector','BackgroundColor','white');
    fprintf('Saved image PDF: %s\n',pdfFile);
    clear guard;
end
end

function close_figure(f)
if isgraphics(f), close(f); end
end
function rr=take(records,im,panel,kind,backend)
rr={};
for k=1:numel(records)
    z=records{k};
    if z.image==im && strcmp(z.panel,panel) && (isempty(kind)||strcmp(z.method.kind,kind)) && ...
            (isempty(backend)||strcmp(z.backend,backend)), rr{end+1}=z; end %#ok<AGROW>
end
end
function rr=setting(records,im,kind,a,rho)
rr=take(records,im,'sensitivity',kind,'');
rr=rr(cellfun(@(z)z.method.alpha==a && z.method.rho==rho,rr));
end
function m=choice(selected,im,kind)
m=[];
for j=1:size(selected,2), v=selected{im,j}; if ~isempty(v)&&strcmp(v.kind,kind), m=v; return; end, end
end
function tf=same_method(a,b)
tf=strcmp(a.kind,b.kind) && a.alpha==b.alpha && a.rho==b.rho && a.tau==b.tau && a.sigma==b.sigma;
end
function a=values(rr,field), a=cellfun(@(z)z.(field),rr); end
function a=finalvals(rr,field), a=cellfun(@(z)z.final.(field),rr); end
function s=count_success(rr), s=sprintf('%d/%d',sum(cellfun(@(z)z.success,rr)),numel(rr)); end
function s=pm(x,d)
if isempty(x), s='--'; return; end
if any(~isfinite(x)), s='n/a'; return; end
if numel(x)==1, sd=0; else, sd=std(x,0); end
s=sprintf('$%.*f \\pm %.*f$',d,mean(x),d,sd);
end
function y=median_all(x)
if isempty(x)||any(isnan(x)), y=NaN; else, y=median(x); end
end
function s=sci(x)
if isnan(x), s='n/a'; elseif isinf(x)
    if x>0, s='$\infty$'; else, s='$-\infty$'; end
elseif x==0, s='$0$'; else
    e=floor(log10(abs(x))); a=x/10^e; s=sprintf('$%.2f\\times10^{%d}$',a,e);
end
end
function s=number(x,d), s=sprintf('%.*g',d,x); end
function s=seed_string(v)
s=strtrim(sprintf('%d, ',v)); if ~isempty(s)&&s(end)==',', s(end)=[]; end
end
function s=plain_name(s)
s=method_name(s); s=strrep(s,'Condat--V{\~u}','Condat-Vu'); s=strrep(s,'--','-');
end
function s=method_name(s)
switch s, case 'MRI',s='MRI-NFB'; case 'RI',s='RI-NFB'; case 'CV',s='Condat--V{\~u}'; case 'CP',s='Chambolle--Pock'; end
end
function s=tex(x)
x=char(x); s='';
for k=1:numel(x)
    switch x(k)
        case '\', a='\textbackslash{}';
        case {'_','%','&','#','$','{','}'}, a=['\' x(k)];
        case '^', a='\textasciicircum{}'; case '~', a='\textasciitilde{}';
        otherwise, a=x(k);
    end
    s=[s a]; %#ok<AGROW>
end
end
function row(fid,parts), fprintf(fid,'%s \\\\\n',strjoin(parts,' & ')); end
function fid=openout(path)
fid=fopen(path,'w'); if fid<0, error('MRI:WriteFailed','Cannot write %s.',path); end
end
