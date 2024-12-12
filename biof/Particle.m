classdef Particle<handle
    %PARTICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius;
        position;
        velocity;
        deposed=0;
    end
    
    methods
        function obj=Particle(radius,position,velocity) %constructor
            if(nargin>0)
                obj.radius=radius;
                obj.position=position;
                obj.velocity=velocity;
            else
                obj.radius=1;
                obj.position=[0,0,10];
                obj.velocity=[0,0,-1];
            end
        end
        function obj=collide(obj) %describe what happens during a collision here
            obj.velocity=zeros(3,1);
            %obj.deposed=deposed;
          
        end
        
        function w=inside2d(obj,point)
            w=0;
            if (norm([obj.position(1)-point(1);obj.position(2)-point(2)])<obj.radius)                
                w=1;
            end
        end
    end
    methods (Static)
        function obj=newInstance(radius,position,velocity)
            obj=Particle(radius,position,velocity);
        end
    end 
end

