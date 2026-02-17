# segment-measure-classify-cells
Macro to quantify the fluorescent signals in the nuclei/cells for all channels and manually classify the cells (pos/neg) based on log-thresholding
 
Brief workflow:
- Nuclei/cells are detected using the StarDist/Cellpose convolutional neural network model for fluorescent nuclei.
- Nuclear ROIs can be filtered on size and eroded to prevent edge effects.
- Background signal in the measurement channel is measured (several options).
  N.B. Intensity values shown in the Results table are already background-subtracted.
- Classifcation can be performed on each channel separately. A histogram of the logarithm of the mean/median intensity or stddev is shown.
  The user sets the threshold by clicking in the histogram image

Input: a folder containing 2D images with at least 2 channels. Multi-series proprietary microscopy format files (e.g. .czi, .lif) are also supported.

Required Fiji update sites:
- StarDist
- CSBDeep
- TensorFlow
- CLIJ
- CLIJ2
- SCF MPI CBG
- PTBIOP

 N.B. This script heavily relies on (mostly) CLIJ2 (https://clij.github.io/), StarDist (https://imagej.net/StarDist) and Cellpose (https://github.com/MouseLand/cellpose).
 If you use this script in a publication, please cite them appropriately.
