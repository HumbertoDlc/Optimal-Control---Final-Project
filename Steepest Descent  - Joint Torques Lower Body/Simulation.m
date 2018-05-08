%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                         OPTIMAL CONTROL                         %%%%%
%%%%%                                                                 %%%%%
%%%%%                                                                 %%%%%
%%%%%                          Final Project                          %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Simulation based on the dyamics obtained by optimization:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% States:
load('x.mat')
load('par.mat')

for k=1:1:10
q(k,:)=interp(x(k,:),5);
end
clear x

qd = zeros(size(q));
    qdd = zeros(size(q));
    xfw = 0;
    nframes = length(q);
    
    % Prepare the new file
    vidObj = VideoWriter('tmp.avi');
    open(vidObj);
    
    for j=1:1:3
    for i = 1:nframes
        % use the rowerdynamics function to get the stick figure
        % coordinates at frame i
        [~,~,~,~,~,~,~,~,s] = rowerdynamics(q(:,i)',qd,qdd,xfw,par);
        clf
        plot(s(:,1),s(:,2),'o-');
        hold on
        
        % draw the seat path
        x = [-1 1];
        y = x*par.a + par.b;
        plot(x,y,'LineWidth',2);
        axis([-1 1 0 1.5]);%'equal');
        
 % Write each frame to the file
        if i==1
            currFrame = getframe;
        else
            % to make sure the frame is always the same size
            currFrame = getframe(gca,[0 0 vidObj.Width vidObj.Height]);
        end
        writeVideo(vidObj,currFrame);
        %pause(0.01)
    end
    end
% Close the file.
 close(vidObj);
     