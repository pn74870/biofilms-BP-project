function [ cov,v,positions,deposed ] = Simulation1( step_size,time,rate )
%Simulation of spherical particle deposition

tic; %start timer
v=0;
persistent particleRadius;
persistent xMax;
persistent xMin;
persistent yMin;
persistent yMax;
persistent N;


positions=zeros(time*rate,3);
deposed=zeros(1,time*rate);

particleRadius=1e-6;
overlaps=zeros(1,10);
xMax=1e-5;
xMin=-1e-5;
yMax=1e-5;
yMin=-1e-5;
zMax=2e-5;
zMin=2e-5;
%velocity=[0;0;-50];

f=gobjects(1,time*rate); %uncomment to draw each frame

N=0; %number of created particles

calcCov=1;


%---------------
T=300; %K
d=1e-6;
eta=1e-5;
mass=1e-15;
gforce=mass*[0,0,-9.8];
%---------------
nOfZSteps=100;
figure();
hold on;
plane1=[xMin,xMax,xMax,xMin];
plane2=[yMin,yMin,yMax,yMax];
plane3=[0,0,0,0];

fill3(plane1,plane2,plane3,'r');
alpha 0.1;
axis equal;
axis([2*xMin,2*xMax,2*yMin,2*yMax,0,zMax]);
cov=zeros(1,time);

interval=1/rate;
nIntervals=time*rate;

for t=1:nIntervals
    
    
        positions(t,:)=[xMin+rand()*(xMax-xMin);yMin+rand()*(yMax-yMin);zMin+rand()*(zMax-zMin)];
        N=N+1;
        %physics simulation
        
        physicsStep(interval);
            
            
%             while (any(particle.velocity))
%                 if(checkForCollision(particle))
%                     particle.collide();
%                     delLastParticle=0;
%                 end
%                 if (particle.position(3)<particle.radius)
%                     particle.collide();
%                     particle.deposed=1;
%                     delLastParticle=0;
%                 end
%                 [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,fluidVel);
%                 particle.position=particle.position+[dx;dy;dz];
%                 %  lastParticle=drawSphere(particle.position);
%                 %  delLastParticle=1;
%             end
            
    %drawSphere(particle.position);
    if (mod(t,rate)==0 && calcCov)
        cov(t/rate)=coverage1();
    end
    fprintf('%.1f percent completed\n',t/nIntervals*100);
end

drawAll();

if (calcCov)
    figure();
    plot(cov);
end
%volume fraction as a function of z
% zTop=max(positions(:,3));
% zStep=zTop/nOfZSteps;
% v=zeros(1,nOfZSteps);
% zVolF=linspace(0,zTop,nOfZSteps);
% for s=1:nOfZSteps
%     v(s)=volumeFraction(zStep*(s-1)+particleRadius,zStep);
% end

% figure();
% plot(zVolF,v);

 fprintf('time %.3f s\n',toc);

    function drawAll()
        for i=1:N
            [x,y,z]=sphere;
            surf(x*particleRadius+positions(i,1),y*particleRadius+positions(i,2),z*particleRadius+positions(i,3));
            
        end
    end
%    function c=coverage()
%        totalArea=(xMax-xMin)*(yMax-yMin);
%        overlapArea=0;
%        nIn=0;
%        inside(1,N)=Particle;
%        for j=1:N
%            within(xMin,xMax,particles(j).position(1));
%            within(yMin,yMax,particles(j).position(2));
%            if(within(xMin,xMax,particles(j).position(1)) && within(yMin,yMax,particles(j).position(2)))
%                nIn=nIn+1;
%                inside(nIn)=particles(nIn);
%               % text(inside(nIn).position(1),inside(nIn).position(2),inside(nIn).position(3)+particleRadius*1.5,'in');
%            end
%        end
%        
%        for i=1:nIn
%           for k=1:nIn
%               if (i~=k)
%                 dist=norm([inside(i).position(1);inside(i).position(2)]-[inside(k).position(1);inside(k).position(2)]);
%                 q=dist/2/particleRadius;
%                 if(q<=1)
%                     overlapArea=overlapArea+particleRadius.^2*(acos(q)-q*sqrt(1-q.^2));
%                 end
%               end
%           end
%        end
%        c=(nIn*pi*particleRadius^2-overlapArea)/totalArea;
%        end
%    
       function w=within(min, max, x)
           w=0;
           if(x>min && x<max)
               w=1;
           end
       end
   
 function shift=periodicBoundary(partIndex)
        sideX=xMax-xMin;
        sideY=yMax-yMin;
        shift=[0,0,0];
        if(xMin>positions(partIndex,1))
            shift(1)=sideX;
        elseif (xMax<positions(partIndex,1))
            shift(1)=-sideX;
        end
        
        if(yMin>positions(partIndex,2))
            shift(2)=sideY;
        elseif (yMax<positions(partIndex,2))
            shift(2)=-sideY;
        end
        
    end
            
   function c=coverage1()
        totalArea=(xMax-xMin+2*particleRadius)*(yMax-yMin+2*particleRadius);
        coveredArea=0;
       % g=zeros(numbOfSteps,numbOfSteps);
        parfor i=1:N
            if (positions(i,3)<=particleRadius)
                coveredArea=coveredArea+pi*particleRadius^2;
            end
        end
        c=coveredArea/totalArea;
   end
            
            
            
    function c=coverage(numbOfSteps)
        totalArea=(xMax-xMin)*(yMax-yMin);
        stepX=(xMax-xMin)/numbOfSteps;
        stepY=(yMax-yMin)/numbOfSteps;
        coveredArea=totalArea;
       % g=zeros(numbOfSteps,numbOfSteps);
        parfor i=1:numbOfSteps
           
            for j=1:numbOfSteps
                free=1;
     
                for par=1:N
                    w=0;
                    point=[xMin+i*stepX,yMin+j*stepY];
                    if (norm([positions(par,1)-point(1);positions(par,2)-point(2)])<particleRadius)
                        w=1;
                    end
                    if(w && deposed(par))
                       free=0;
                       break;
                     %  g(mod(i,numbOfSteps)+1,floor(i/numbOfSteps)+1)=1;
                    end
                end
                if(free)
                    coveredArea=coveredArea-stepX*stepY;
                end
            end
        end
        c=coveredArea/totalArea;
%         figure();
%         hold on;
%         for x=1:numbOfSteps
%             for y=1:numbOfSteps
%                 if (g(x,y)==1)
%                     fill3([(x-1)*stepX,x*stepX,x*stepX,(x-1)*stepX],[(y-1)*stepY,(y-1)*stepY,y*stepY,y*stepY],[0,0,0,0],'b');
%                 end
%             end
%         end
%         hold off;
    end

function fig=drawSphere(p)
    
        if(~isempty(f))
            try
                delete(f(p));
            catch
            end
        end
        [x,y,z]=sphere;
        fig=surf(x*particleRadius+positions(p,1),y*particleRadius+positions(p,2),z*particleRadius+positions(p,3));
        drawnow;
        
    end


    function v=volumeFraction(z,thickness)
        totalVol=(xMax-xMin)*(yMax-yMin)*thickness;
        occVol=0;
        for i=1:N
            if(deposed(i) && positions(i,3)<z+thickness+particleRadius  && positions(i,3)>z-particleRadius) %within the box
                h1=thickness-positions(i,3)+z;
                h2=positions(i,3)-z;
                capVol1=sphCapVolume(particleRadius,h1);
                capVol2=sphCapVolume(particleRadius,h2);
                occVol=occVol+4*pi/3*particleRadius^3-capVol1-capVol2;
            end
        end
        v=occVol/totalVol;
        %-------------------------------------
        %nested function to calculate the volume of spherical cap
        function v=sphCapVolume(R,h)
            if(h<R)
                v=pi*(2/3*R^3-R^2*h+h^3/3);
            else
                v=0;
            end
        end
        %-----------------------------------------
    end

    function overlap=checkCol(i)
        overlap=0;
                    for j=1:N
                        %collision detection
                        if(norm(positions(j,:)-positions(i,:))<2*particleRadius && i~=j) %overlapping spheres
                            overlap=1;
                            break;
                        end
                    end
    end
    function physicsStep(time)
        while (time>0)
            for i=1:N
                if(~deposed(i))
                    
                    [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,gforce);
                    move=[dx,dy,dz];
                    move=move+periodicBoundary(i);
                    positions(i,:)=positions(i,:)+move;
                   
                        while (checkCol(i))
                            positions(i,:)=positions(i,:)-move./10;
                            
                        
                        end
                    
                    %                             dist=positions(j,:)-positions(i,:);
                    %                             if(minD>dist)
                    %                                 minD=dist;
                    %                             end
                    %-------------------------------------------------
                    %   if(deposed(j)) %collision with a fixed particle
                    %      resolveCollision(i,j);
                    %                             if(positions(j,3)<=particleRadius)
                    %                                 resolveCollision(i,j);
                    %                             else %collision with a moving particle
                    %                                 dCenters=positions(j,:)-positions(i,:);
                    %                                 positions(i,:)=positions(i,:)-0.5*(2*particleRadius-norm(dCenters))*normr(dCenters);
                    %                                 positions(j,:)=positions(j,:)+0.5*(2*particleRadius-norm(dCenters))*normr(dCenters);
                    %
                    %                             end
                    %                          %-------------------------------------------------
                    
                end
                
            
            
            
            
            if(positions(i,3)<particleRadius) %collision with the surface
                positions(i,3)=particleRadius;
                deposed(i)=1;
                
            end
            f(i)=drawSphere(i); %uncomment to draw each frame
            
            end
            
            time=time-step_size;
        end
end
        
    function resolveCollision(movingPartIndex,stationaryPartIndex)
       
        dCenters=positions(stationaryPartIndex,:)-positions(movingPartIndex,:);
        positions(movingPartIndex,:)=positions(movingPartIndex,:)-(2*particleRadius-norm(dCenters))*normr(dCenters);
        %deposed(movingPartIndex)=1;
    end
   
end


