This new repository replaces the microtrack repository.

We are trying to clean out the old and unused routines, get ready for building the fascicle evaluation structures, and to write simple demonstration functions and scripts that will be the foundation for more complex experiments.

FP, BW

OLD COmments

We are developing the simpler form of microtrack here.

Nov. 1, 2011
The original and complex code from Tony is in the sherbondy-old directory.

Here we started to develop scripts for the NIH grant that used the same ideas but in simpler coding form.

The first pass on the new microtrack calculations is in the directory arcuate. In there we analyzed STT fibers and STT plus contrack fibers to see how many we would be able to delete by culling (dtiCullFiber) followed by the microtrack analysis in which we build a matrix to predict the diffusion signal and then remove the fibers with very small weights (near zero).

The resulting fibers ended up predicting the diffusion data reasonably well.  We are now moving on to other computational experiments.

The namespace for microtrack scripts is s_mct<ScriptName> or mct<FunctionName>

BW, FP, JY
