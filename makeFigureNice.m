function [] = makeFigureNice( fig_size, fig_pos )
% simple function to make figure look nicer.
% 
% Inputs:
%     fig_size        - OPTIONAL
%                       Choose dimensions of figure [x size, y size]
%     fig_pos         - OPTIONAL
%                       Choose position of figure [x pos, y pos]

if nargin == 0
   fig_size = [600, 400];
end
if nargin < 2
   fig_pos = [100, 100];
end

% Size the figure
fig = gcf;
set(fig,'Position',[fig_pos(1), fig_pos(2), fig_size(1), fig_size(2)]);

% Set the fonts
set(gca,'FontSize',16);
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'XLabel'),'FontSize',16);
set(get(gca,'Title'),'FontSize',16);

% Make the lines thicker
hPlot = get(get(fig,'Children'),'Children');

try     % won't execute if there are no lines to thicken
    if iscell(hPlot)
        set(hPlot{1},'LineWidth',2.0);
        set(hPlot{2},'LineWidth',2.0);
    else
        set(hPlot,'LineWidth',2.0);
    end
catch
    % do nothing
end

% turn the grid on
grid on;

% make tick marks thicker
ax = gca;
ax.TickLength       = 3*ax.TickLength;
ax.LineWidth        = 1.5;
ax.GridLineStyle    = '--';

% % Copy to clipboard
% print -dbitmap -r200;

end