addpath ../mytoolbox/geo_tools/

r_dir=strcat(' ../../data_refine/'); %read directory
s_dir=strcat(' ../../data_refine/'); %save directory

periodicBc=[true,false];
shearlineOdeSolverOptions = odeset('relTol',1e-4);
strainlineOdeSolverOptions = odeset('relTol',1e-2);
forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];
strainLcsColor= [1,0,0];

for iii=1:20
    load_ftle=strcat('load ',r_dir,'ftle_field',num2str(iii),'.mat');
    load_ow=strcat('load ',r_dir,'ow_field',num2str(iii),'.mat');
    
    eval(load_ftle);
    run.flow.resolution=resolution;
    run.flow.domain=domain;

    cgEigenvalue=cgStrainD;
    cgEigenvector=cgStrainV;
    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    cgEigenvector1 = reshape(cgEigenvector(:,1:2),[fliplr(resolution),2]);  
    
    imagesc(domain(1,:),domain(2,:),dle); hold on
    colormap(flipud(gray))
    drawnow
    
    rcgStrain=permute(reshape(cgStrain,2,2,prod(resolution)),[3 4 1 2]);
    poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
    
    eval(load_ow);
    pSe=poincareSection(4).endPosition;
    poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
    poincareSection(1).endPosition=pSe;

    hPoincareSection = arrayfun(@(input)plot(input.endPosition(:,1),input.endPosition(:,2)),poincareSection);
    set(hPoincareSection,'color',shearLcsColor)
    set(hPoincareSection,'LineStyle','--')
    set(hPoincareSection,'marker','o')
    set(hPoincareSection,'MarkerFaceColor',shearLcsColor)
    set(hPoincareSection,'MarkerEdgeColor','w')
    drawnow

    [poincareSection.numPoints] = deal(45);
    nPoincareSection = numel(poincareSection);
    for i = 1:nPoincareSection
        rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
        poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
    end

%    lamb_vec=linspace(max(cgStrainD(:,1)),min(cgStrainD(:,2)),30);%<----should the number of layers be variable?
    lamb_vec=linspace(1,1.25,60);
    
    for lam_i=1:length(lamb_vec)
        saver = strcat('save ',s_dir,'geo_',num2str(iii),'_lambda',num2str(lam_i+60),'.mat closedOrbits dle lamb_vec domain resolution');
        disp('Detect elliptic LCS ...')
        
        lambda = lamb_vec(lam_i);
        [shearline.etaPos,shearline.etaNeg] = lambda_line(cgEigenvector,cgEigenvalue,lambda);
        closedOrbits = poincare_closed_orbit_multi(run.flow,shearline,poincareSection,'odeSolverOptions',shearlineOdeSolverOptions,'showGraph',false,'periodicBc',periodicBc);

        hClosedOrbitsEtaPos = arrayfun(@(i)plot(closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
        set(hClosedOrbitsEtaPos,'color',shearLcsColor)
        set(hClosedOrbitsEtaPos,'linewidth',2)
        hClosedOrbitsEtaNeg = arrayfun(@(i)plot(closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
        set(hClosedOrbitsEtaNeg,'color',shearLcsColor)
        set(hClosedOrbitsEtaNeg,'linewidth',2)
        for j = 1:nPoincareSection
            hOrbitsPos = cellfun(@(position)plot(position(:,1),position(:,2)),closedOrbits{j}{1});
            set(hOrbitsPos,'color',strainLcsColor)
            hOrbitsNeg = cellfun(@(position)plot(position(:,1),position(:,2)),closedOrbits{j}{2});
            set(hOrbitsNeg,'color',strainLcsColor)
        end
        drawnow
        eval(saver)
    end
end