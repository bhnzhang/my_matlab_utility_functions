function [] = f_add_path_with_subfolders( top_level_folder )
% Adds path with subfolders
%
% if is a git repository, removes .git folders from the path
% actually i haven't figured out how to do that yet
%
% Inputs:
%   top_level_folder
%       type: string
%       desc: top level folder

addpath( genpath( top_level_folder ) );

end

