function [main] = shadederrbar(x,y,sem,varargin)
%Name: shadederrbar
%Created By: Pratik Mistry
%Created On: 27 July 2018
%Updated: 31 July 2018 -- Turned the main plot into a output variable.
%
%[main] = shadederrbar(x,y,sem,color)
%
%Description: This code creates a plot of the main graph with transparent
%shaded error bars when it is provided with the x vector, y vector,
%standard error mean, and a color
%Input: x-vector, y-vector, SEM, RGB color value (**all row vectors**)
%Output: A graph
%
%Usage: 
%Note: If you would like to create multiple plots, type figure before calling
%this function.
%
% figure
% shadederrbar(time,velocity,SEM,[0 0 1])
% figure
% shadederrbar(days,temperature,SEM,[1 0 1])
%
%
%Note: If you like to create plots on the same graph, just call the
%function again.
%
%figure
%shadederrbar(time,velocity,SEM,[0 1 1])
%shadederrbar(time,acceleration,SEM,[0 1 0])

    if nargin == 4
        color = varargin{1};
    else
        color = [0 0 0]; %Default plot color is black 
    end

    if(size(x,1)~=1); x = x'; end
    if(size(y,1)~=1); y = y'; end
    if(size(sem,1)~=1); sem = sem'; end
    
    if (ischar(color)==1) %Checks to see if the color inputted is a numerical RGB value
        color=char2rgb(color);
    end
    patchcolor=color+(1-color)*.8; %Creates the patch color
    
    yerru=y+sem;
    yerrl=y-sem;
    
    xpatch=[x,fliplr(x)]; %Creates x axis for the path
    ypatch=[yerru,fliplr(yerrl)]; %Creates y axis for patch
    
    hold on
    % plot(x,yerru,'-','Color',patchcolor); %Plots upper error line
    % plot(x,yerrl,'-','Color',patchcolor); %Plots lower error line
    fill(xpatch,ypatch,patchcolor,'FaceAlpha',0.5,'EdgeAlpha',0,'EdgeColor',patchcolor,'LineStyle','none'); %Creates the patch
        main = plot(x,y,'-','Color',color,'HandleVisibility','off'); %Plots main data   
    % hold off
    
end
