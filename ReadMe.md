# Linear Fascicle Evaluation (LiFE)

Standard tractography can use diffusion measurements from a living brain to generate a large collection of candidate white-matter fascicles; the connectome. Linear Fascicle Evaluation (LiFE) takes any connectome and uses a forward modelling approach to predict diffusion measurements in the same brain. LiFE predicts the measured diffusion signal using the orientation of the fascicles present in a connectome. LiFE uses the difference between the measured and predicted diffusion signals to measure prediction error. The connectome model prediction error is used to compute two metrics to evaluate the evidence supporting properties of the connectome. One metric -the strength of evidence - compares the mean prediction error between alternative hypotheses. The second metric - the earth movers distance - compares full distributions of prediction error. These metrics can be used for: 1. Comparing tractography algorithms 2. Evaluating the quality of tractography solutions for individual brains or group of brains and 3. Testing hypotheses about white-matter tracts and connections.

### License.
#### Copyright (2013-2014), Franco Pestilli, pestillifranco@gmail.com

### Application.
* Evaluate the evidence supporting white-matter connectomes generated using [magnetic resonance diffusion imaging](http://en.wikipedia.org/wiki/Diffusion_MRI) and [computational tractography ](http://en.wikipedia.org/wiki/Tractography).

* Perform statistical inference on white-matter connectomes: Compare white-matter connectomes, show the evidence for white-matter tracts and connections between brain areas.

### Installation.
1. Download [LiFE](https://github.com/vistalab/life).
2. Download [vistasoft](https://github.com/vistalab/vistasoft).
3. Download [MBA](https://github.com/francopestilli/mba).
4. [Start MatLab](http://www.mathworks.com/help/matlab/startup-and-shutdown.html).
5. Add LiFE and vistasoft to the [matlab search path](http://www.mathworks.com/help/matlab/ref/addpath.html).

### Dependencies.
1. [MatLab](http://www.mathworks.com/products/matlab/).
2. [vistasoft](https://github.com/vistalab/vistasoft).
3. [Matlab Brain Anatomy (MBA)](https://github.com/francopestilli/mba).

### [Documentation.](http://vistalab.github.io/life/doc/)

### Getting started.
1. [Demo MatLab File](http://vistalab.github.io/life/doc/Pestilli_etal_manuscript/life_demo.html)
2. [Demo Results](http://vistalab.github.io/life/html/life_demo.html)
3. [Demo data](http://purl.stanford.edu/cs392kv3054)

### [Stable code release.](https://github.com/vistalab/life/releases/tag/v0.2)

### How to cite LiFE.
Pestilli F., Yeatman J.D., Rokem A., Kay K.N., Wandell B.A. Linear fascicle evaluation (LIFE) of white matter connectomes. Poster presentation at the Organization for Human Brain Mapping Annual Meeting, Seattle, WA, June 2013.
