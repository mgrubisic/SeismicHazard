function[shakefield]=buildshake(model,scenario,opt)

Nscenarios = size(scenario,1);

shakefield(1:Nscenarios,1)=struct('label',[],'datasource',[],'type',[],'mechanism',[],'surface',[],'vertices',[],'thickness',[],'gptr',[],'geom',[],'gmpe',[],'mscl',[],'rupt',[],'mulogIM',[],'L',[],'Z',[],'gpsRA',[]);

for i=1:Nscenarios
    ptr1 = str2double(scenario{i,1}); % model pointer
    slib = {model(ptr1).source.label};
    [~,ptr2] = intersect(slib,scenario{i,2});
    source = model(ptr1).source(ptr2);
    
    % updates source
    source.mscl.M   = scenario{i,6};
    source.mscl.dPm = 1;
    
    rot   = source.geom.rot;
    pmean = source.geom.pmean;
    hypmr = [scenario{i,7},scenario{i,8},0];
    hypmr0 = bsxfun(@minus,source.geom.hypm,pmean)*rot;
    hypmr0(:,3)=[];
    d2     = sum(bsxfun(@minus,hypmr0,hypmr(1:2)).^2,2);
    [~,ind]=min(d2);
    hypm  = bsxfun(@plus,hypmr*rot',pmean);
    
    source.geom.aream  = 0;
    source.geom.hypm   = hypm;
    source.geom.normal = source.geom.normal(ind,:);
    source.mulogIM =[];
    source.L   =[];
    source.Z   =[];
    
    % Computes Rupture Area coordinates (for graphical Purpose Only)
    M = scenario{i,6};
    X = scenario{i,7};
    Y = scenario{i,8};
    RAtype = scenario{i,5};
    rupt = source.rupt;
    switch RAtype
        case 'circular'
            Amax = rupRelation(M,0,Inf,rupt.RA{1},rupt.RA{2});
            r    = sqrt(Amax/pi);
            t    = (0:pi/30:2*pi)';
            xyz   = [X+r*cos(t),Y+r*sin(t),t*0];
            xyzRA =  bsxfun(@plus,xyz*rot',pmean);
        case 'rectangular'
            Amax = rupRelation(M,0,Inf,rupt.RA{1},rupt.RA{2});
            aratio = rupt.aratio;
            L    = sqrt(Amax*aratio);
            W    = sqrt(Amax/aratio);
            dX   = [-1 1 1 -1]'*L/2;
            dY   = [-1 -1 1 1]'*W/2;
            xyz   = [X+dX,Y+dY,dX*0];
            xyzRA =  bsxfun(@plus,xyz*rot',pmean);
        case 'adaptive'
    end
    
    gpsRA = xyz2gps(xyzRA,opt.ellipsoid);
    source.gpsRA  = gpsRA;
    shakefield(i)=source;
end

