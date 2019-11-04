function [] = makeFigureNice( OPTS )
% simple function to make figure look nicer.
% 
% Inputs:
%     OPTS
%       type: struct
%       desc: OPTIONAL options
%             'grid_line_style': 
%                   any acceptable matlab grid line style,
%                   such as '--' for dashed or '-' for solid
%             'fig_size':
%                       Choose dimensions of figure [x size, y size]
%             'fig_pos':
%                       Choose position of figure [x pos, y pos]

% DEFAULTS
if nargin == 0
   OPTS.fig_size        = [800, 600];
   OPTS.fig_pos         = [100, 100];
   OPTS.grid_line_style = '--';  % default to dashed lines
end
if nargin < 2
   
end
if nargin < 3
   OPTS.grid_line_style = '--';  % default to dashed lines
end

% Size the figure
fig = gcf;
set(fig,'Position',[OPTS.fig_pos(1), OPTS.fig_pos(2), OPTS.fig_size(1), OPTS.fig_size(2)]);

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
ax.GridLineStyle    = OPTS.grid_line_style;

% % Copy to clipboard
% print -dbitmap -r200;

end