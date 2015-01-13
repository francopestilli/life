# Linear Fascicle Evaluation (LiFE)

Standard tractography can use diffusion measurements from a living brain to generate a large collection of candidate white-matter fascicles; the connectome. Linear Fascicle Evaluation (LiFE) takes any connectome and uses a forward modelling approach to predict diffusion measurements in the same brain. LiFE predicts the measured diffusion signal using the orientation of the fascicles present in a connectome. LiFE uses the difference between the measured and predicted diffusion signals to measure prediction error. The connectome model prediction error is used to compute two metrics to evaluate the evidence supporting properties of the connectome. One metric -the strength of evidence - compares the mean prediction error between alternative hypotheses. The second metric - the earth movers distance - compares full distributions of prediction error. These metrics can be used for: 1. Comparing tractography algorithms 2. Evaluating the quality of tractography solutions for individual brains or group of brains and 3. Testing hypotheses about white-matter tracts and connections.

## Application.
* Evaluate the evidence supporting white-matter connectomes generated using [magnetic resonance diffusion imaging](http://en.wikipedia.org/wiki/Diffusion_MRI) and [computational tractography ](http://en.wikipedia.org/wiki/Tractography).

* Perform statistical inference on white-matter connectomes: Compare white-matter connectomes, show the evidence for white-matter tracts and connections between brain areas.

## License.
#### Copyright (2013-2015), [Franco Pestilli](http://francopestilli.com/), pestillifranco@gmail.com

## [Documentation](http://francopestilli.github.io/life/doc/).

## [Stable code release](https://github.com/vistalab/life/releases/tag/v0.2).

## How to cite LiFE.
[Pestilli, Franco, Jason D. Yeatman, Ariel Rokem, Kendrick N. Kay, and Brian A. Wandell. Evaluation and statistical inference for human connectomes. Nature methods 11, no. 10 (2014): 1058-1063.](http://www.nature.com/nmeth/journal/v11/n10/abs/nmeth.3098.html)

## LiFE in [Python](https://www.python.org/)[:Dipy](http://nipy.org/dipy/)
The LiFE algorithm has been recently implemented in Python by [Ariel Rokem](http://arokem.org/) and is now available as part of the Dipy software: [LiFE@Dipy](http://nipy.org/dipy/examples_built/linear_fascicle_evaluation.html#example-linear-fascicle-evaluation).

## Installation.
1. Download [LiFE](https://github.com/francopestilli/life).
2. [Start MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).
3. Add LiFE to the [matlab search path](http://www.mathworks.com/help/matlab/ref/addpath.html).

## Dependencies.
* [MatLab](http://www.mathworks.com/products/matlab/).
* [vistasoft](https://github.com/vistalab/vistasoft).
* [Matlab Brain Anatomy (MBA)](https://github.com/francopestilli/mba).

## Getting started.
Learn about LiFE by using [life_demo.m](http://francopestilli.github.io/life/doc/scripts/life_demo.html) in [MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).

### 1. [Download LiFE](https://github.com/francopestilli/life).
* Download the LiFE repository from the TAR/ZIP files linked [here](https://github.com/francopestilli/life/archive/v0.2.zip).
* UNZIP/UNTAR the file.
* Add the life folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/life/folder/'))
```

### 2. [Download vistasoft](https://github.com/vistalab/vistasoft).
* Download the VISTASOFT repository from the TAR/ZIP files linked [here](https://github.com/vistalab/vistasoft/archive/master.zip).
* UNZIP/UNTAR the file.
* Add the VISTASOFT folder to your matlab search path. To do so in the MatLab prompt type: 
```
   >> addpath(genpath('/my/path/to/the/VISTASOFT/folder/'))
```

### 3. [Download LiFE Data Demo](http://purl.stanford.edu/cs392kv3054).
* Download the LiFE demo data set from the repository [here](https://stacks.stanford.edu/file/druid:cs392kv3054/life_demo_data.tar.gz).
* UNZIP/UNTAR the file.
* Add the unzipped/untarred Data folder to your matlab search path. To do so in the MatLab prompt type:
```
   >> addpath(genpath('/my/path/to/the/life_data_demo/folder/'))
```

### 4. [Read the life_demo documentation](http://vistalab.github.io/life/doc/scripts/life_demo.html).
Read the description of the calculations in the documentation inside the file, life_demo.m by typing the following in the matlab prompt: 
```
  >>  edit life_demo
```

### 5. [Run the life_demo code](https://github.com/francopestilli/life/blob/master/scripts/life_demo.m).
This final step will run the life_demo code. The code will perform the operations described [here](http://vistalab.github.io/life/html/life_demo.html). 
```
  >>  life_demo
```
life_demo.m runs in about 30 minutes on a modern Intel processor with 8GB of RAM. This code has been tested with MatLab 2012b on Ubuntu 12.10 and Mac OSX 10.9.


