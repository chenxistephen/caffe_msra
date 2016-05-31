function build(debugbuild, varargin)

clear mex

vtr = getenv('VisionToolsRoot');
if (isempty(vtr))
    error('VisionToolsRoot environment variable not set');
end

compflags = ['COMPFLAGS="$COMPFLAGS"'];
debugflags = 'DEBUGFLAGS="/Zi"';
linkflags = ['LINKFLAGS="$LINKFLAGS /LTCG /NODEFAULTLIB:"libcmt.lib""'];

cnfg = 'release';
miscflags = []; 
if (exist('debugbuild', 'var') && debugbuild ~= 0)
    miscflags = ['-g OPTIMFLAGS="/Od"'];
    cnfg = 'debug'
    compflags = [compflags ' -D_DEBUG'];
    linkflags = ...
        ['LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:"libcmtd.lib""'];
end

incdir = [vtr '\inc'];
libdir = [vtr '\lib'];
dsrcdir = [vtr '\src\detection'];
if (strcmp(mexext, 'mexw64'))
    libs = ['-lcommon' cnfg 'x64 -lnumerics' cnfg 'x64' ...
        ' -lcore' cnfg 'x64 -ldetection' cnfg 'x64'];
else
    libs = ['-lcommon' cnfg 'win32 -lnumerics' cnfg 'win32' ...
        ' -lcore' cnfg 'win32 -ldetection' cnfg 'win32'];
end

outpath = [fileparts(mfilename('fullpath')) '\'];

% build channels
output = [outpath 'vtchannels' '.' mexext];
mexstr = ['mex' ' ' dsrcdir '\channelsmex.cpp' ' ' ...
	dsrcdir '\vtmatlab.cpp' ' ' ...
    '-I' incdir ' ' '-L' libdir ' ' ...
    libs ' ' compflags ' ' linkflags ' ' debugflags ' ' ...
    '-output' ' ' output ' ' '-v' ' ' miscflags]
eval(mexstr)

% build detector
output = [outpath 'vtdetector' '.' mexext];
mexstr = ['mex' ' ' dsrcdir '\detectormex.cpp' ' ' ...
	dsrcdir '\vtmatlab.cpp' ' ' ...
    '-I' incdir ' ' '-L' libdir ' ' ...
    libs ' ' compflags ' ' linkflags ' ' debugflags ' ' ...
    '-output' ' ' output ' ' '-v' ' ' miscflags]
eval(mexstr)

% build helpers
output = [outpath 'vthelpers' '.' mexext];
mexstr = ['mex' ' ' dsrcdir '\helpersmex.cpp' ' ' ...
	dsrcdir '\vtmatlab.cpp' ' ' ...
    '-I' incdir ' ' '-L' libdir ' ' ...
    libs ' ' compflags ' ' linkflags ' ' debugflags ' ' ...
    '-output' ' ' output ' ' '-v' ' ' miscflags]
eval(mexstr)
