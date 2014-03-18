%{
Run this file to compile the mex files if you need to
(or if the binaries work for your computer, then you can skip this)
This file also contains two simple test problems

I tested this on linux, so it might not work out-of-the-box on Windows or Mac OS,
but I don't foresee any major difficulties in getting it to work. 

The main issue with Windows is finding a Fortran compiler. If you don't
have one, then you can try the C version (which I generated with f2c),
which has worked sometimes for some users. 

Stephen Becker, Feb 14, 2012  srbecker@alumni.caltech.edu
Updated         May  3, 2012  adding the f2c version for Windows userss
Updated         Jun 24, 2012  adding a lot of trouble-shooting support mainly for 64-bit linux


See also lbfgsb.m and lbfgsb_wrapper.c

%}

if ispc
    disp('Try the compile_mex_forWindows.m script, especially if you don''t have a Fortran compiler');
end

%% -- Download lbfgsb v. 3.0 from Jorge Nocedal's website -- 
% From his website: Condition for Use: This
% software is freely available, but we expect that all publications describing  work
% using this software , or all commercial products using it, quote at least one of
% the references given below. This software is released under the BSD License

url ='http://users.eecs.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz';
srcDir = untar(url,'./');

%% -- Compile --

% run mex -setup if you haven't already...

% With mex, you can 'define' certain symbols by using the -D command,
% and 'undefine' them with the -U command, e.g.
%   -DDEBUG will define the "DEBUG" symbol, and -UDEBUG will undefine it
%   (if defined, it gives more verbose output when you run the program )

if ispc
    mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
        Lbfgsb.3.0/lbfgsb.f Lbfgsb.3.0/linpack.f Lbfgsb.3.0/timer.f ...
        -lm -lmwblas  COMPFLAGS="$COMPFLAGS -O3"
    
    % You might also have better luck by first compiling the .f files
    %   then run mex -setup to choose a C compiler, and then compile
    %   and link the .c and .obj files
else
    % For Linux and Mac.
    % I haven't tested this in Darwin/Mac OS, but it might work
    
    % Note: using g95 as the Fortran compiler doesn't seem to work
    %   I had much better luck with gfortran (version 4.3)
    %   The FC="gfortran" part of the command will force Matlab to use gfortran
    %   (you can also edit the settings in prefdir/mexopts.sh
    
    mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
        Lbfgsb.3.0/lbfgsb.f Lbfgsb.3.0/linpack.f Lbfgsb.3.0/timer.f ...
        -lm -lblas  CFLAGS="\$CFLAGS -O3" FC="gfortran"

end

%% Troubleshooting


% BLAS problems
%   BLAS is the library that does all matrix multiplications
%   Matlab comes with their own library, called libmwblas, and it
%   is usually very fast (it incorporates Intel or AMD's libraries)
%   and it can usually take advantage of multicore processors. So
%   try that.  The compilation command is:
%       -lmwblas
%   You might also want to try -D_BLAS64_
%
%   Use the BLAS that comes with your system by specifying
%       -lblas
%
%   Unfortunately, not all BLAS is compatible, since there are 32-bit and
%   64-bit flavors. So try several BLAS if it doesn't work at first.
%   You can also modify the source code by hand (if you want to modify
%   C code instead of Fortran code, then see the f2c'ed version
%   which I describe in compile_mex_forWindows.m )
%
%   Matlab has a 32-bit compatible version of BLAS too. To use this,
%   the command is
%       -lmwblascompat32 blas32header.h 

% Problems specific to Windows
%   When compiling both Fortran and C programs, the names of functions
%   are changed differently depending on what compiler you have.
%   So the main Fortran function, setulb, can be known to the linker
%   as SETULB, or setulb_, or SETULB_
%
%   The flags to pland around with
%       -DNOUNDERSCORE or -UNOUNDERSCORE
%       -DUPPERCASE_FORTRAN or -UUPPERCASE_FORTRAN
%       -D_WIN32 or -U_WIN32

if ~ispc && strfind(computer, '64' )
    % for 64-bit on linux and mac
    
    % You can try this. It is using -lmwblas instead of -lblas. However,
    % if you compiled with -lblas and then the program crashed when you
    % tried to run it, it will probably also crash when you use -lmwblas.

    %mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
        %Lbfgsb.3.0/lbfgsb.f Lbfgsb.3.0/linpack.f Lbfgsb.3.0/timer.f ...
        %-lm -lmwblas  CFLAGS="\$CFLAGS -O3" FC="gfortran"

    % Instead, I recommend using Matlab's 32-bit compatible BLAS.

    % First, we run a little script that will call 'sed' and rename
    %   the blas functions to their equivalent 32-bit version
    % (e.g. ddot is renamed to ddot32)
    [status,result] = system('./convertTo32.sh');
    if status
        error('script failed; try running script by hand perhaps?');
    end

    % Now compile, using the 32-bit version
    mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
        Lbfgsb.3.0/lbfgsb32.f Lbfgsb.3.0/linpack32.f ... % using new files
        Lbfgsb.3.0/timer.f ...
        -lm -lmwblascompat32 ...
        CFLAGS="\$CFLAGS -O3" FC="gfortran"

    %   The issue is because recent versions of BLAS
    %   have some 64-bit integer types. You can verify if this is the case
    %   by examining both:
    % fullfile(matlabroot,'extern','include','blas.h')
    %   and 
    % fullfile(matlabroot,'extern','include','blascompat32.h')
    %   and see if they are different. For example, on my 64-bit laptop
    %   with R2010a, blascompat32.h defines the first integer argument to ddot32
    %   as "int", whereas blas.h defines the first integer argument to ddot
    %   as "ptrdiff_t" (which is usually the same as size_t).

end

%% Troubleshooting part 2: trying the f2c version
%   Avoids fortran completely and uses the C version that was automatically
%   (and then modified by hand) converted from the Fortran.
%   I have had limited success with this.

TRY_F2C = false;
% TRY_F2C = true;

if TRY_F2C
    
    if ~exist('lbfgsb.c','file')
        % Download it from my website
        % (as of June 2012, the .c files should be included in the main filecentral package )
        url     = 'http://ugcs.caltech.edu/~srbecker/code/Lbfgsb_f2c.3.0.tgz';
        srcDir  = untar(url,'./');
        disp('Done downloading files');
    end
    
    % adding this flag might help you compile:
    %   -UUSE_TIMER
    if ispc
        disp('See compile_mex_forWindows.m');
    else
        % Assuming that you already have the f2c library installed
        %   You can also use -lf2c and Matlab will use the dynamic library
        %   (e.g. libf2c.so for linux, libf2c.dylib for mac)
        %   but then you might have problems at runtime with complaints
        %   about a missing reference to MAIN__. So I suggest using the static
        %   library.
        
        % Try default BLAS:
        if false
          mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
            lbfgsb.c linpack.c timer.c -I/opt/matlab/extern/include/...
            -lm -lblas -U_BLAS64_ CFLAGS="\$CFLAGS -O3" ... 
            /usr/lib/libf2c.a
        
        else
          % Try Mathworks' better BLAS:
          if false
            %   (but it may have 64-bit bindings which can be a pain)
            mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
                lbfgsb.c linpack.c timer.c ...
                -lm -lmwblas -D_BLAS64_ -U_BLAS32_ CFLAGS="\$CFLAGS -O3" ...
                /usr/lib/libf2c.a 
          else
            % If this fails, then you can try -lblascompat32.
            % Since the code is now in C, it's a bit easier, and
            % we don't have to run the bash script to edit the files
            % Instead, just pass in -U_BLAS64_ and -D_BLAS32_
            mex -v lbfgsb_wrapper.c -largeArrayDims  -UDEBUG ...
                lbfgsb.c linpack.c timer.c ...
                -lm -lmwblascompat32 -U_BLAS64_ -D_BLAS32_ CFLAGS="\$CFLAGS -O3" ...
                /usr/lib/libf2c.a
          end
        end
    end
end
%% test the new function
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
