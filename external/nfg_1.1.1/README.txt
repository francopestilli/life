----- Numerical Fibre Generator v1.1 (29/4/09) ------

This document contains practical information regarding the use of the Numerical Fibre Generator.  For a more detailed 
description of the technique and theory behind it please refer to the related paper, which at the time of writing this
document was in press.  Also read the descriptions provided in the default parameters files for more understanding on how 
the tools can be modified.

Thomas G. Close, Jacques-Donald Tournier, Fernando Calamante, Leigh A. Johnston, Iven Mareels and Alan Connelly, A 
software tool to generate simulated white matter structures for the assessment of fibre-tracking algorithms,
 NeuroImage (2009), doi:10.1016/j.neuroimage.2009.03.07/7
   

Contents:

 1. License Info (for full license see 'GNU_GENERAL_PUBLIC_LICENSE.txt')
 
 2. Citation Info (if you use this software please note the recommended citation format below)
 
 3. Overview

 4. Data Structure

 5. File Formats

	- Data
	- Parameters
	- DW Gradient Directions/b-values
	- MR Images

 6. Tool Descriptions

	- rand_init
	- optimise
	- subdiv
	- trim
	- resample
	- mri_sim
	- noisify
	- draw_rois

 7. Compiling

 8. Running
 
 9. Displaying

	- display_strands
	- display_subvoxels
	- display_colour_key
	- display_stats
	- display_ext_stats
	- display_rois

----------------
1. License Info
----------------

   Created by Tom Close on 25/06/08.
   Copyright 2008 Tom Close.
   Distributed under the GNU General Public Licence.
 
   'Numerical Fibre Generator' is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   'Numerical Fibre Generator' is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with 'Numerical Fibre Generator'.  If not, see <http://www.gnu.org/licenses/>


----------------
2. Citation Info
----------------

If you use my software in a publication could you please include a line such as the following to cite my work: 

"The synthetic white matter and the associated DW-MR images were created using the 'Numerical Fibre Generator' software package 
(T. G. Close, Brain Research Institute, Melbourne, Australia, http://wwww.brain.org.au/software/) [1]." (see citation below)

At the time of writing this document the following paper was in press, 

[1] Thomas G. Close, Jacques-Donald Tournier, Fernando Calamante, Leigh A. Johnston, Iven Mareels and Alan Connelly, A 
software tool to generate simulated white matter structures for the assessment of fibre-tracking algorithms,
 NeuroImage (2009), doi:10.1016/j.neuroimage.2009.03.07/7

-------------
3. Overview
-------------

The numerical fibre generator is a collection of seven 'command-line' tools: 'rand_init', 'optimise', 'subdiv', 'trim', 'resample',
 'mri_sim' and 'noisify'.  It is used to simulate the structural characteristics of human white matter and the diffusion-weighted MR 
images that would arise from them.

Each tool performs one stage of processing.  They are usually used in the sequence above (with the exception that 'optimise' 
and 'trim' are usually used twice).  However each tool can be used independently on appropriately formatted datasets.  

------------------
4. Data Structure
------------------

The white matter structures are represented as collections of continuous 'strands'.  These strands are defined as continuous 
paths with nominal cross-sectional radii.  These paths are the linear interpolation of a sequence of control points.  

The start and end points are fixed (and are hence considered separately), while the rest of the control points are free to move 
during the optimisation (see 'optimise' tool description below).  In order to define the orientation a strand at its start and 
end points, a 'pre point' and 'post point' are prepended and appended to the strand respectively.  The line joining the pre 
(post) and start (end) points defines the starting (ending) orientation of the strand.

If included fixed (spherical) isotropic regions are defined by centre coordinate, a radius about this centre, a modelled diffusivity (mm/s),
a signal intensity relative to the b=0 intensity of the strands, and a relative repulsion weighting.  The relative repulsion weighting determines the
relative weighting in the cost function of overlap between the isotropic region and surrounding strands.  This weighting corresponds to the length
of strand (of equal radius) required to produce the same overlap repulsion.


------------------
5. File Formats
------------------

--- Strands ---

The file format used to hold the strands is quite straightforward.  A 'strand collection' (the output at each stage) is stored in 
a single directory containing a seperate text file for each strand.  These files are named in the following format:


	'strand_' + (strand index) + '-' + (strand bundle index) + '-r' + (strand cross-sectional radius) + '.txt'


The pre, start, control, end and post points are stored in the 'strand file', in that order.  Each point is stored on a seperate 
line with their three coordinates (x,y,z) seperated by single spaces.  

The strand radii, coordinate positions, voxel size, ...  are only relative to one another.  Thus they can represent any scale
that is convenient (eg. mm, cm, ...).


--- Isotropic Regions ---

All isotropic regions are stored within a single file located within the strand collection directory called 'isotropic_regions.txt'.
The parameters for each region are stored separate lines in the following format

x_coordinate	y_coordinate	z_coordinate	radius	diffusivity	relative_b=0_intensity	relative_repulsion_weighting


--- Parameters ---

Parameters are stored in text files, with each parameter stored on a separate line.  The format for each line is simply:

key value

The parameters 'key' (name) must be at the very beginning of the line, followed by a SINGLE space then the parameter's value
corresponding value.  Any lines that do not start with a recognised parameter are ignored. 

WARNING: Lines that start with a recognised key and then are not followed by a valid value (or a valid value not after a SINGLE 
space) will not be checked!  


--- DW Gradient Directions / b-values ---

DW gradient directions and their corresponding b values are also stored in text files.  Each gradient direction/b-value pair is 
stored on a separate line.  The gradient direction is stored in cartesian coordinates (x,y,z) and is normalised once it is 
read in. The format for each pair is:

X Y Z b-value

Each value is separated by a SINGLE space and the X coordinate must be right at the beginning of the line.


--- MR Images ---

The MR images simulated by NFG are stored in ANALYZE format, with the voxels ordered by increasing coordinates in the x, y then z 
dimensions. 

NB: While both image simulation tools ('mri_sim' and 'noisify') accept an 'output_format' parameter, the only recognised option 
at this stage is 'analyze'.  If you want to add a different image format then please use this parameter to signal it.   


--------------------------------------
6. Tool Descriptions
--------------------------------------

Descriptions of what each tool does is given below.  Each tool has a number of parameters that control aspects of its function.  
Default parameters can be found in appropriately named text files (eg. 'rand_init_param.txt', 'optimise_param.txt' ... ) in the 
directory './default_parameters' along with descriptions as to which aspect each parameter is controlling.

NB: For a more detailed understanding on how each tool works I suggest reading the parameter file comments as almost all parameters
are tuneable and the comments descriptions on the aspects they are controlling. 


--- rand_init ---

Usage:  arg[1]= output directory location
	arg[2]= parameters file location

'rand_init' randomly initialises a collection of strands such that their start and end points lie on the surface of a sphere.  
It randomly assigns each strand a radius and ensures that no strands overlap at the sphere surface.  The strand radii are drawn 
from a uniform distribution between a predefined upper and lower bound.  However if, given this radius, the start or end points 
of the strand overlap with a previously generated strand the new strand is rejected.  This means that radii closer to the lower 
bound are increasingly likely as the surface of the sphere is filled, thus skewing the distribution of strand radii towards the 
lower bound.  

Additionally, to ensure that the centre of the sphere is sufficiently filled, longer strands are preferred via a rejection 
sampling scheme where new strands are only accepted if 

	dot_product(start_point, end_point) < X

Where X is a uniformly distributed random variable between -1 and 1.  



--- optimise ---

Usage:  arg[1]= input directory location
	arg[2]= output directory location
	arg[3]= parameters file location

'optimise' is the heart of the toolbox, it is the tool that ensures that there is minimal overlap between strands while 
maintaining reasonable curvature and length.  

It optimises the coordinates of the control points, while keeping the start and end points fixed, w.r.t. a cost function 
(source code located at './c_source/optimise/cost_function.c') that penalises overlap between strands, strand length and strand 
curvature.  Depending on the complexity of the initial configuration the optimisation algorithm will usually converge to an 
acceptable configuration within 100 iterations.  Due to the difficulty in defining suitable stopping criteria that are 
invariant with different numbers/sizes of strands, the optimisation algorithm is left to run until it can no longer reduce the 
cost function. 

After each iteration the new 'cost', its calculated gradient, the step size taken and so forth are displayed to screen, in 
addition to a '.' being printed after every cost function call.  When the optimisation can proceed no further it will stop and 
save the new strand configuration.



--- subdiv ---

Usage:  arg[1]= input directory location
	arg[2]= output directory location
	arg[3]= parameters file location

'subdiv' divides all the strands in the collection into thinner strands of the same radius, specified by an input parameter.  
The new strands will be hexagonally packed across the cross-section of the original strands.




--- trim ---

Usage: 	arg[1]= input directory location
	arg[2]= output directory location
	arg[3]= parameters file location

'trim' trims all control points that lie outside of a specified radius (given as an input parameter).  Strands that enter the 
sphere more than once are split into multiple strands. The start (end) points of the strands are recalculated to be the point 
where the strand enters (leaves) the sphere.   


--- resample ---

Usage: 	arg[1]= input directory location
	arg[2]= output directory location
	arg[3]= parameters file location

'resample' resamples the control points of all the strands in the collection to a new constant interval (given as an input parameter).
This is primarily used to remove bunching and stretching introduced in the subdivision stage. However, it could be used to increase
the sampling frequency before a final fine-tuning optimisation stage.  

The subdivision can also sometimes leave folds in the outer strands, where a successive control point is placed behind (relative to the
parent strand direction) the previous control point.  This occurs when there are moderate kinks in the parent strand and the bundle 
radius is large.  In this case, the moderate kinks in the parent strand are amplified in the distal child strands, which can cause 
the distal strands to temporarily double back on themselves.  These unwanted kinks are checked for by a parameter, 
'double_back_angle_threshold', which determines if the angular deviation at any given control point is due to this effect (i.e. if the 
angular deviation exceeds this threshold).  Once a strand is determined to double back on itself, all successive control points are 
discarded until the angular deviation between a succeeding control point and the control point where the strand was deemed to 
double back on itself falls below a second threshold, 'forward_angle_threshold'.


--- mri_sim ---

Usage:	arg[1]= input directory location
	arg[2]= output file location (extensions for header and data files will be appended)
	arg[3]= imaging gradient directions (normalised cartesian coordinates) and b values 
	arg[4]= parameters file location

'mri_sim' simulates the DW-MR images that would arise from imaging the simulated white matter structure.  Orientations of the 
strands are first plotted onto a fine subvoxel grid to retain partial volume effects.  A tensor model of diffusion is assumed
 at each subvoxel element and the corresponding DW signal is calculated given a specified b-value.  These signals are then 
summed into the coarser voxel grid for the final DW-MR images (ie retaining the partial volume contribution to the signal within each voxel).  

In addition to the images a statistical summary of the strand collection will be written to a file 'stats.txt' in a 
subdirectory 'stats' of the strand collection directory.  Note, images will be in Analyze 7.5 format.


--- noisify ---

Usage:	arg[1]= input directory location
	arg[2]= output directory location
	arg[3]= parameters file location

'noisify' adds rician noise to the images with a standard deviation specified by an input parameter 'noise_level'.  


NB: After running each tool the parameters passed to it and previous tools used on the data set are copied into a subdirectory 
of the strand collection directory called 'parameters'


--- draw_rois ---

Usage:	arg[1]= input directory location
	arg[2]= output directory location
	arg[3]= parameters file location

'draw_rois' draws ROIs (regions of interest) around the start and end of each bundle. These ROIs can then be used as seed and target regions 
for conventional tracking algorithms.  By default all the masks will be saved in one combined signed integer integer, however by setting the 
'save_combined_mask' parameter to 0 the masks will be saved in seperate unsigned char images.  In the combined image, masks are differentiated by
different integers derived from their bundle index and whether they are the start or end ROI by the formula 

				2 * (bundle_index),      		if start ROI
				2 * (bundle_index) + 1,  		if end ROI

where blank voxels are denoted by -1. In the seperate images, the masks are denoted by 1 and the blank voxels 0. The seperate masks will be named as 
follows

	'mask-(bundle_index)-(0 if start or 1 if end).img'.
	
So the mask of the end of of the 15th bundle will be saved in a file

	'mask-15-1.img'
	
and the start of the 8th bundle will be saved in the file

	'mask-08-0.img'
	
Note, that not all bundles will make it to the final structure as some will be completely trimmed. These bundles will not have a mask generated for them. 
The masks will be in Analyze 7.5 format.


--------------
7. Compiling
--------------


All tools are known to compile with gcc v4.0.  Before compiling you will need to install the GNU Scientific Library v1.9 (or 
greater).  It can be found at 'www.gnu.org/software/gsl/' (downloaded at ftp://ftp.gnu.org/gnu/gsl/).

There is a makefile in the base directory.  You may need to alter the 'GSL_PATH' and 'CC' variables in the makefile to suit your 
installation of GSL and gcc respectively.  Simply type 'make' to compile all the tools (NB: remember to type 'make clean' before making any 
subsequent attempts).   

After compiling the executables can be found in the './bin' directory.  In order to use the 'run' bash script you should add it to your path. 

--- Windows ---

NFG has been successfully compiled on both 'MinGW/MSYS' and 'Cygwin'. If you don't have either then I would recommend MinGW/MSYS (http://www.mingw.org/)
as it is simpler/smaller.  They can be downloaded from "http://sourceforge.net/projects/mingw/".  Download the MinGW automated installer and the
current release of MSYS (but install MinGW first). 

NB: Install these first before installing GNU Scientific Library v1.9.



--------------
8. Running
--------------

A bash script called './run' will call each of the tools in an appropriate sequence to simulate a white matter structure and two sets of noisy and no noise 
DW-MRI images with 20 and 60 diffusion gradient directions respectively. 
 
 Usuage:
 
	arg[1]= (optional) Parameter directory
 
		By default each tool will use the parameters stored in separate, appropriately named files in the './parameters/default' directory.  
		Parameter files used can be changed by passing the name of a different parameters directory as the first argument to the run script.
		Only the parameter files that have been altered are required in the new directory, if the script doesn't find the parameters file 
		for a given tool it will use the default.
	
	arg[2]= (optional) Gradient directions file location
 
		By default the 60 directions, 3000 b value file './gradient_directions/encoding_60dir.txt' will be used.  Use this argument to 
		specify a different gradient directions/b value file.  
	
	arg[3]= (optional) Base output directory

		By default the generated strand collections will be stored in subdirectories of the './generated_collections' directory.  Use this
		argument to specify a different base output directory. New collections will be created in subdirectories named accordingly to the 
		date and time they are created (see argument 4).
	
	arg[4]= (optional) Output directory	
	
		By default the output from each run is stored in a directory named after the datetime the strand collection was initialised 
		(in 'yymmddHHMMSS' format). The strand collection at each processing stage will be stored in appropriately named subdirectories with the 
		final sets of images being stored in '9_noisify-20dir' and '9_noisify-60dir'.  


NB: For good performance with default settings (performance is approximately inversly proportional to (number of strands) * (sample density)) 
NFG requires at least ~ 1Gb of ram.


--- Windows ---

In order to run the 'run' script you will need to install 'MSYS' or 'Cygwin' (http://www.mingw.org/), which is a Bourne shell environment for Windows.
If you don't have either then I would recommend MSYS (http://www.mingw.org/) as it is simpler/smaller.  It can be downloaded from "http://sourceforge.net/
projects/mingw/".  Download the current release of MSYS. 


--------------
9. Displaying
--------------

Also included is a MATLAB script used to display the generated strands.

--- display_strands ---

Usage:	arg[1]= Either the directory name of the stored collection or a cell structure containing the loaded strands from a 
			   previous call (return variable 1).
		arg[2]= (optional) The number of divisions around the diameter of the strand. The default is 6, which corresponds to 
			   a hexagonal cross-section.
		arg[3]= (optional) A Nx3 matrix containing the colours in which the different bundles will be displayed.
		arg[4]= (optional) The indices of the strands to be displayed.

Output:	return[1]= A 4xN cell structure containing the strand data (to be used as arg[1] in subsequent calls of the function)
		return[2]= A Nx3 matrix containing the colours of the strands displayed (to be used as arg[3] in subsequent calls of 
				 the function to keep colours consistent)

Displays a 3D-rendered image of the strands in the collection.  


--- display_subvoxels ---

Usage:	arg[1]= The directory name of output of a 'mri_sim' call (needs the 'save_subvoxels' parameter to be set to 1 when 
			   calling 'mri_sim').
		arg[2]= (optional) The out of plane slice dimension
		arg[3]= (optional) The slice number.

Displays a RGB direction-encoded fibre orientation plot at each of the subvoxels used to get the final voxel intensities.  If 
no value for slice number is supplied it will iterate through all the slides (press any key to move to next slice).  Note that 
the 'save_subvoxels' parameter needs to be set to 1 when calling 'mri_sim' and that if the resolution is reasonably high (eg. 
300 x 300 x 300) then the file size can be quite large (> 600MB). 

--- display_colour_key --

Usuage: arg[1] = An N x 3 'colour' matrix (as returned in the second output of 'display_strands')

Displays an image that can be used as a legend for the different bundle colours displayed in 'display_strands'


--- display_stats --- 

Usuage: arg[1] = The directory name of output of a 'mri_sim' call (needs the 'save_subvoxels' parameter to be set to 1 when 
			   calling 'mri_sim').

Simply reads the 'stats/stats.txt' and displays it in the matlab command window


--- display_ext_stats ---

Usuage: arg[1] = The directory name of output of a 'mri_sim' call (needs the 'save_ext_segment_stats' parameter to be set to 1 when 
			   calling 'mri_sim').
	vararg[2, ..., 5] = (optional) Either 'segment_length', 'segment_angle', 'segment_overlap_fraction' or 'segment_radius_curv'
				 specifying the extended characteristic to display.  If no arguments supplied all will be displayed.

Displays histograms of the specified characteristics. 
  

--- display_rois ---

Usage:	arg[1]= Directory name of the stored collection 
		arg[2]= (optional) A Nx3 matrix containing the colours in which the different bundles will be displayed.
		arg[3]= (optional) The out of plane slice dimension
		arg[4]= (optional) The slice number (if left blank or set to zero then all slices will displayed in succession).

Displays the rois as drawn by the tool 'draw_rois'.  Each loaded roi is coloured according to its bundle index. Can uses the same 'colour_key'
as 'display_strands' so it is possible to compare their figures. Note that by altering the default parameters of 'draw_rois' it is possible to produce 
masks of complete bundles.

Requires the 'save_combined_mask' parameter of 'draw_rois' to be set to 1.
