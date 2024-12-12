function [ cov ] = Simulation( step_size,time,rate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

particleRadius=1e-6;
xMax=1e-5;
xMin=-1e-5;
yMax=1e-5;
yMin=-1e-5;
zMax=5e-5;
zMin=5e-5;
velocity=[0;0;-50];
particles(1,500)=Particle;
N=0; %number of created particles


delLastParticle=0; %used for animation

%---------------
T=30;
d=1e-6;
eta=1e-5;
fluidVel=[0,0,-1e-4];
%---------------

figure();
hold on;
plane1=[xMin,xMax,xMax,xMin];
plane2=[yMin,yMin,yMax,yMax];
plane3=[0,0,0,0];
fill3(plane1,plane2,plane3,'r');
axis equal;
axis([2*xMin,2*xMax,2*yMin,2*yMax,0,zMax]);
cov=zeros(1,time);
for t=1:time
    for p=1:rate
        particle=Particle.newInstance(particleRadius,[xMin+rand()*(xMax-xMin);yMin+rand()*(yMax-yMin);zMin+rand()*(zMax-zMin)],velocity);
        particles((t-1)*rate+p)=particle;
        N=N+1;
        %physics simulation
        
        while (any(particle.velocity))
            if(checkForCollision(particle))        
                particle.collide();
                delLastParticle=0;
            end
            if (particle.position(3)<particle.radius)
                particle.collide();
                particle.deposed=1;
                delLastParticle=0;
            end
            [dx,dy,dz]=BrownianMotion(step_size,T,eta,d,fluidVel);
            particle.position=particle.position+[dx;dy;dz];
          %  lastParticle=drawSphere(particle.position);
          %  delLastParticle=1;
        end
        drawSphere(particle.position);
        
    end
    cov(t)=coverage(100);
    
end
figure();
plot(cov);

function c=checkForCollision(particle)
    c=0;
    for j=1:N
        if (norm(particles(j).position-particle.position)<2*particleRadius && ~isequal(particles(j),particle) && particles(j).deposed)
            c=1;
            particle.deposed=1;
            
        end
    end
end

    function fig=drawSphere(pos)
        if(delLastParticle)
            delete(lastParticle);
        end
        [x,y,z]=sphere;
        fig=surf(x*particleRadius+pos(1),y*particleRadius+pos(2),z*particleRadius+pos(3));
        drawnow;
        
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
%        function w=within(min, max, x)
%            w=0;
%            if(x>min && x<max)
%                w=1;
%            end
%        end
    function c=coverage(numbOfsteps)
        totalArea=(xMax-xMin)*(yMax-yMin);
        stepX=(xMax-xMin)/numbOfsteps;
        stepY=(yMax-yMin)/numbOfsteps;
        coveredArea=totalArea;
        g=ones(numbOfsteps,numbOfsteps);
        for i=1:numbOfsteps
            for j=1:numbOfsteps
                free=1;
                for par=1:N
                    if(particles(par).inside2d([i*stepX;j*stepY]))
                       free=0;
                       g(i,j)=0;
                    end
                end
                if(free)
                    coveredArea=coveredArea-stepX*stepY;
                end
            end
        end
        c=coveredArea/totalArea;
        figure();
        title('jk');
        hold on;
        for i=1:numbOfsteps
            for j=1:numbOfsteps
                if (g(i,j)==1)
                    ra = stepX/2;
                    theta = 0:0.1:6.29;
                    x = i*stepX + ra*cos(theta);
                    y = j*stepY + ra*sin(theta);
                    
                    fill(x,y,'b')
                    plot(x,y,'b')
                end
            end
        end
    end
end


