%{

Stephen Becker, April 17, 2012  srbecker@alumni.caltech.edu

Extra installation script for Windows, assuming you do not have
a Fortran compiler.

The lbfgsb package is in Fortran, and compiling Fortran on Windows
is not fun nor easy. So I have used the "f2c" utility to convert
the Fortran files to C files. The only problem is that "f2c" relies
on its own libraries, which are a little tricky to get working.
You can download precompiled f2c libraries (e.g. "libf2c.lib") from the
web, but these have been compiled with specific compilers. If your Matlab
compiler doesn't match those versions, then you might have issues.

So, this script compiles the f2c library from scratch. The lbfgsb C files
might not need the entire library, so my approach is to include the entire
library and then exclude files that cause problems. So far, it seems
to work.

If compilation fails, add "-v" flags to the mex scripts and you will
get more verbose and useful output.

%}

if ~ispc
    disp('Run the other install script; on linux, installation should be simple');
end

%% Step -1: download the precompiled binary from my website and see if it works
% Download the file and skip directly to the tests at the end
% Make sure you have saved anything important, because there's a chance
% that running this file will crash Matlab.

if strcmpi(computer,'pcwin')
    url     = 'http://ugcs.caltech.edu/~srbecker/code/lbfgsb_wrapper.mexw32';
    untar(url,'./');
elseif strcmpi(computer,'pcwin64')
    % Sorry, I don't have a precompiled binary.
    % If you manage to compile it, please email it to me
    % and I will host it on my website. stephen.beckr@gmail.com
    disp('Sorry');
    %url     = 'http://ugcs.caltech.edu/~srbecker/code/lbfgsb_wrapper.mexw64';
    %untar(url,'./');
else
    disp('If you''re not running windows, then don''t run this sript');
end

%% Step 0: run mex -setup if you haven't already.
% This tells Matlab which compiler to use

%   32-bit Matlab on Windows comes with the lcc-win32 compiler, which is
%   what I tested this with.

%   64-bit Matlab doesn't come with a compiler, but you can get a free
%   Microsoft C compiler; see
%   http://www.mathworks.com/support/compilers/R2011b/win64.html

% mex -setup

%% Step 1, old. Windows users should skip this step unless they have a Fortran compiler

% Download lbfgsb v. 3.0 from Jorge Nocedal's website -- 
% From his website: Condition for Use: This
% software is freely available, but we expect that all publications describing  work
% using this software , or all commercial products using it, quote at least one of
% the references given below. This software is released under the BSD License

% url ='http://users.eecs.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz';
% srcDir = untar(url,'./');

%% Step 1, new and improved (for Windows users):
%   Download the f2c version of lbfgsb v3.0
%   I have converted this from Fortran to C using f2c, and made a few small
%   changes by hand. It depends on the f2c libraries, which we will
%   deal with shortly.

url     = 'http://ugcs.caltech.edu/~srbecker/code/Lbfgsb_f2c.3.0.tgz';
srcDir  = untar(url,'./');
% This zip file contains lbfgsb.c, linpack.c, and timer.c

%% Step 2: Download f2c library files
BASE=pwd;
url         ='http://www.netlib.org/f2c/libf2c.zip';
srcLocation = fullfile(BASE,'f2c_src');
srcDir      = unzip(url, fullfile(BASE,'f2c_src') );
%% Step 2b: Rename a few of the .h0 files to .h
cd(srcLocation)
for fileC = {'f2c','sysdep1','signal1'}
    file = fileC{1};
if ~exist([file,'.h'],'file')
    if exist([file,'.h0'],'file')
        success = movefile([file,'.h0'],[file,'.h']);
        if ~success, disp('Please manually copy the files'); end
    else
        error('Cannot find the files; are you in the right directory?');
    end
else
    fprintf('Found file %s.h; no need to do anything\n', file );
end
end
%% Step 3: Try to compile all the f2c files using Matlab's compiler
cd(srcLocation)
list = dir('*.c');
% Make the .obj files, but don't link yet
for j = 1:length(list)
    src_file = list(j).name;
    % Some C files don't compile with lcc, but they don't appear to be essential
    switch src_file
        case {'ftell64_.c','pow_qq.c','qbitbits.c','qbitshft.c', ...
            'signbit.c','uninit.c' }
            fprintf('Skipping %s\n', src_file);
        otherwise
            if ~exist( [src_file(1:end-1),'obj'],'file' )
                fprintf('Compiling %s\n', src_file );
                mex('-largeArrayDims','-I./','-c','-DMSDOS','-DNOUNDERSCORE','-DNO_ISATTY',src_file);
            end
    end
end
objList = dir('*.obj');
%% Step 4: Compile the LBFGSB files, my LBFGSB wrapper file, and link
cellList =  {objList.name};
cellList_withDir =  {};  % a list of the f2c object files
for j = 1:length(cellList)
    obj = cellList{j};
    % exclude certain files, like main.obj
    switch obj
        case {'main.obj','arithchk.obj','derf_.obj','derfc_.obj',...
                'getarg_.obj','s_paus.obj','erf_.obj','erfc_.obj','iargc_.obj'}
            fprintf('Skipping %s\n', obj );
        otherwise
            cellList_withDir{end+1} = ['f2c_src',filesep,obj];
    end
end
cd(BASE);

mex('-v','-largeArrayDims','-If2c_src',...
    'lbfgsb.c','linpack.c','timer.c', ...
    '-output','lbfgsb_wrapper', ...
    'lbfgsb_wrapper.c','-lmwblas', cellList_withDir{:} );

disp('Done compiling');

%% And now we are done. Run some tests below to verify it works

%% Simple test

disp('=== lbfgsb "driver1" test problem, 2D === ');
% Here's the test problem included with lbfgsb called 'driver1'

n   = 25;

l   = ones(n,1); u = l;
odd = 1:2:n;
even= 2:2:n;
l(odd) = 1.0;
u(odd) = 1.0e2;
l(even)= -1.0e2;
u(even)=  1.0e2;


opts    = struct( 'x0', 3*ones(n,1) );
opts.printEvery     = 1;
opts.m  = 5;

[x,f,info] = lbfgsb( @driver1, l, u, opts );

% The true objective value is 0.
if abs(f) < 1e-8
    disp('Success!');
else
    disp('Something didn''t work right :-(  ');
end

% the structure info.err contains the objective function (1st column)
%   and norm(gradient,Inf) (2nd column)

semilogy( info.err(:,1)-f,'o-' ); xlabel('iteration'); ylabel('error in objective function');

%% another test function, the 2D Rosenbrock function
disp('=== Rosenbrock test function, 2D === ');
n = 2;

fxy = @(x,y) 100*( y-x.^2).^2  +  (1-x ).^2 ;
f   = @(x)   fxy( x(1,:), x(2,:) );
gxy = @(x,y) [100*(4*x.^3-4*x.*y)+2*x-2; 100*(2*y-2*x.^2)];
g   = @(x)   gxy( x(1,:), x(2,:) );

% There are no constraints
l   = -inf(n,1);
u   = inf(n,1);

opts    = struct( 'x0', [-1.9;2] );
opts.printEvery     = 1;
opts.m  = 5;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = [1;1];
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1e3;

% The {f,g} is another way to call it
[x,f,info] = lbfgsb( {f,g} , l, u, opts );

if abs(f) < 1e-8
    disp('Success!');
else
    disp('Something didn''t work right :-(  ');
end

% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
