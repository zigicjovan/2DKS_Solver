function [outFile, vwo] = gif2avi(file, destination, varargin)
% gif2avi(file) converts gif image 'file' to an avi video file and saves
% it to the current directory with the same file name. 'file' is a
% character vector or string scalar containing a file name or full path
% to a gif file and must include the .gif extension.
%
% gif2avi(file,destination) is used to specify the path and/or file name
% of the output file in one of the following forms:
%   * 'path\file.ext' 
%   * 'path' - default file name is same as input file name
%   * 'file.ext' - default path is current directory: cd()
%   * '.ext' - default file name and path are used.
% The file extension ext can be .avi or .mp4 but see 'profile' options.
% Destination can be empty ([],'') etc to use default value.  
%
% gif2avi(file,[],'nLoops', n) records the gif loop n-times where n is 
% a positive integer.  Default = 1;
%
% gif2avi(file,[],'profile', str) specifies which profile to use in 
% VideoWriter. See doc('VideoWriter') for options. Profile must be 
% compatible with the output file format (avi or mp4). Default values are
% 'Indexed AVI' for avi output and 'MPEG-4' for mp4 output.
%
% gif2avi(file,[],'FrameRate', fps) sets the frame rate of the output file  
% in frames per second. See VideoWriter documentation for help. 
%
% gif2avi(file,[],'Colormap', m) sets the colormap of the output file. See
% VideoWriter documentation for help. 
%
% gif2avi(file,[],'Quality', s) sets the video qualtiy of the output file 
% between [0,100] where higher values result in larger higher quality files.
% See VideoWriter documentation for help. 
%
% outFile = gif2avi(__) returns the full path to the avi file which can be 
% used in implay(outFile).
% [~, vwo] = gif2avi(__) returns the VideoWriter object.
%
% When complete, a message will be printed to the command window
% displaying the avi file name and path.  For PC users, the message will
% contain a link that will open the directory that stores the avi file.
%
% If the gif file has varying colormaps for each frame (uncommon), the 
% colormapping of the avi file may not be correct due to a limitation in
% writeVideo.
%
% Large GIF images may require a lot of processing time.
%
% EXAMPLES
%   gif2avi('myGifImage.gif')
%   gif2avi('myGifImage.gif, 'myAviImage.avi')
%   gif2avi('myGifImage.gif, '.mp4')
%   gif2avi('myGifImage.gif, 'C:\Users\name\Documents\MATLAB')
%   gif2avi('myGifImage.gif, 'C:\Users\name\Documents\MATLAB\myAviImage.avi')
%   gif2avi('myGifImage.gif, [], 'profile', 'Uncompressed AVI')
%   gif2avi('myGifImage.gif, [], 'FrameRate', 10, 'nLoops', 10)
%   outFile = gif2avi('myGifImage.gif');  implay(outFile)
%
% Source: <a href = "https://www.mathworks.com/matlabcentral/fileexchange/76198-gif2avi">gif2avi()</a>
% Copyright (c) 2020, <a href = "https://www.mathworks.com/matlabcentral/profile/authors/3753776-adam-danz/activities">Adam Danz</a>
% All rights reserved
% gif2avi() was inspired by <a href = "https://www.mathworks.com/matlabcentral/answers/422599">Walter Roberson's answer</a> in the Answers Forum.

% Version history
% 1.0.0     200526  first uploaded to FEX.
% 1.0.1     200526  Updated source in help section to FEX link; added copyright at end.
%                   Added output and now using input file name as default output filename.
%                   Assert that checks extension only accepts avi. Profile input added.
% 1.0.2     200526  Corrected error caused by duplicate perods in file extension. Added
%                   assert to check that destination exists. 
% 1.1.0     200526  fps input added; vwo output added; documentation updated.
% 2.0.0     200527  Partially incompatible with previous versions.   All inputs to VideoWriter
%                   are now name-val pairs replacing inputs 3&4 from prev vs + nLoops.  Now allows 
%                   for mp4 file output.  Input validation changed & additional validation. 

%% Input parser & validation
% Validate 'file' input
narginchk(1,inf)
nargoutchk(0,2)
validateattributes(file,{'char','string'},{'scalartext'},mfilename,'file',1)
% Check that input file exists & is .gif
assert(exist(file,'file')==2, 'The file input could not be found: %s.',file)
[~, inputFilename, inputExt] = fileparts(file);
assert(strcmpi(inputExt,'.gif'), 'The input file must have the .gif file extension; ''%s'' not accepted.', inputExt)

% Validate 'destination' input
if nargin < 2 || isempty(destination)
    destination = cd(); 
end

% Validate and parse optional name-value params
p = inputParser();
p.FunctionName = mfilename;
addParameter(p, 'nloops', 1, @(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'}));
addParameter(p, 'profile', [], @(x)validateattributes(x,{'char','string'},{'scalartext'})); % more validation later.
addParameter(p, 'FrameRate', [], @(x)validateattributes(x,{'numeric'},{'scalar','positive'})); 
addParameter(p, 'Colormap', [], @(x)validateattributes(x,{'numeric'},{'2d','ncols',3})); 
addParameter(p, 'Quality', [], @(x)validateattributes(x,{'numeric'},{'scalar','integer','>=',1,'<=',100})); 
parse(p,varargin{:})

% Fill in missing path|file|ext to destination
[filepath, filename, ext] = fileparts(destination);
missingParts = double([isempty(filepath), isempty(filename), isempty(ext)]); 
switch true
    case isequal(missingParts, [0 0 1]) % only path provided
        filename = inputFilename;
        filepath = destination; 
        if strcmpi(p.Results.profile, 'MPEG-4')
            ext = '.mp4';
        else
            ext = '.avi';
        end
    case isequal(missingParts, [1 0 0]) % only filename.ext provided
        filepath = cd(); 
    case isequal(missingParts, [1 1 0]) % only .ext provided
        filepath = cd(); 
        filename = inputFilename;
    case isequal(missingParts, [1 0 1]) % path and filename indistinguishable, missing extension
        error('Destination can either be a full path to a file, just the path, or a filename with extension.')
end
outFile = fullfile(filepath, [filename,ext]);
% Check that destination dir exists and ext is OK
assert(any(strcmpi(ext,{'.avi','.mp4'})),'Output file extension ''%s'' not accepted. Extension must be .avi or .mp4.', ext)
assert(exist(filepath,'dir')==7, 'Destination directory does not exist: %s', filepath)

% Set default profile
if any(strcmpi(p.UsingDefaults,'profile')) && strcmpi(ext,'.avi')
    profile = 'Indexed AVI';  % Using default profile with avi file
elseif any(strcmpi(p.UsingDefaults,'profile')) && strcmpi(ext,'.mp4')
    profile = 'MPEG-4';       % Using default profile with mp4 file
else
    profile = p.Results.profile; % Using user-defined profile
end

% Confirm output file and profile are compatible
if strcmpi(ext,'.avi')
    acceptedprofiles = {'Motion JPEG AVI', 'Uncompressed AVI', 'Indexed AVI', 'Grayscale AVI'};
elseif strcmpi(ext,'.mp4')
    acceptedprofiles = {'MPEG-4'}; 
end
assert(any(strcmpi(profile, acceptedprofiles)), ...
    'A %s file is not compatible with a ''%s'' profile. Compatible profiles are [%s].',...
    ext, profile, strjoin(acceptedprofiles,', '))

%% Convert gif to avi
% Read in GIF file
try
    [gifImage, map] = imread(file, 'Frames', 'all');
catch ME
    if strcmpi(ME.identifier,'MATLAB:TooManyInputs')
        ME = addCause(ME, MException(ME.identifier,'FILE may not be a gif file.'));
    end
    rethrow(ME)
end
assert(ndims(gifImage)==4 && size(gifImage,3)==1, ...
    'Unexpected imread output size [%s].', num2str(size(gifImage)))

% Create output file
vwo = VideoWriter(outFile, profile);
% Optional VideoWriter name-val pairs
if ~isempty(p.Results.FrameRate)
    vwo.FrameRate = p.Results.FrameRate;
end
if ~isempty(p.Results.Colormap)
    vwo.Colormap = p.Results.Colormap;
end
if ~isempty(p.Results.Quality)
    vwo.Quality = p.Results.Quality;
end
open(vwo)

% Loop through frames and write to movie file, nloop-times
for n = 1:p.Results.nloops
    for i = 1:size(gifImage,4)
        F = im2frame(gifImage(:,:,1,i),map);
        if strcmpi(profile, 'Grayscale AVI')
            F.colormap = []; % required for grayscale avi.
        end
        writeVideo(vwo,F);
    end
end
close(vwo)

% Display link to file in command window.
if ispc
    disp([mfilename,'.m is complete. <a href="matlab: winopen(''',filepath,''') ">Open directory for ',[filename,ext],'</a>', '.'])
else
    fprintf('%s.m is complete. %s%s saved in %s.\n', mfilename, filename, ext, filepath)
end
 
%% Copyright (c) 2020, Adam Danz
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution
% * Neither the name of nor the names of its
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.