clear all
close all

% layers.info file
layerfile = 'C:\Users\bz\Google Drive\research\popovic group\tapeouts\2019_05_45RFSOI\layers.info';

f = fopen(layerfile, 'r');
% f = fopen('new_layers.txt', 'r');
A = textscan(f, '%s %s %s %d %d %d','HeaderLines',2);
fclose(f);

GL1_name = A{1};
layer_name = A{2};
purp = A{3};
layer_no = A{4};
layer_datatype = A{5};
layer_CSD = A{6}; % WTF is this?

s = {};
s{end+1} = '<?xml version="1.0" encoding="utf-8"?>';
s{end+1} = '<layer-properties>';


for i=1:length(GL1_name)

    if strcmp(purp{i}, 'drawing'), purp{i} = 'draw'; end
    if strcmp(purp{i}, 'exclude'), purp{i} = 'excl'; end

    s{end+1} = ' <properties>';
    s{end+1} = sprintf('  <name>%s %s</name>', layer_name{i}, purp{i});
    s{end+1} = sprintf('  <source>%d/%d</source>', layer_no(i), layer_datatype(i));
    s{end+1} = ' </properties>';
end

s{end+1} = ' <name/>';
s{end+1} = '</layer-properties>';

CRLF = [char(13) char(10)];
TXT = '';

for i = 1:length(s)
    TXT = [TXT, s{i}, CRLF];
end

fid = fopen('new_layers.lyp','w','n');
fwrite(fid, TXT);
fclose(fid);

