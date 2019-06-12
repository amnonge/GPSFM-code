# GPSfM: Global Projective SFM Using Algebraic Constraints
on Multi-View Fundamental Matrices (CVPR 2019)

Yoni Kasten* Amnon Geifman* Meirav Galun Ronen Basri (*equal contribution)


MATLAB Mex Version 1.0 (2019-06-06)
Copyright (C) Yoni Kasten and Amnon Geifman, Weizmann Institute, 2019.
Licensed for noncommercial research use only.


## Background

The code retrievs camera matrices and 3D points given a set of pairwise fundamental matrices.

For more information see:

[[arXiv]](https://arxiv.org/pdf/1812.00426.pdf)
Please cite these paper if you use this code in an academic publication.
```
@InProceedings{Kasten_2019_CVPR,
author = {Kasten, Yoni and Geifman, Amnon and Galun, Meirav and Basri, Ronen},
title = {GPSfM: Global Projective SFM Using Algebraic Constraints on Multi-View Fundamental Matrices},
booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
month = {June},
year = {2019}
}
```

## Installation


mex/C++ code:
In order to use the code it is necessary to compile the mex functions.
We supply compiled versions for Windows and Linux.
For compilation, enter the folder GPSFM\compileMex and run compileMexLinux or compileMexWindows for compiling the Ceres solver for the final bundle adjustment.




## Use

We supply end-to-end code for the pipeline described in the paper.
Important function:
```
ProcessScript.m- Script to generate our data fromat from the format of http://www.maths.lth.se/matematiklth/personal/calle/dataset/dataset.html. We already supply the 25 datasets from the paper in our format in the folder: "DataSet Proj"  If you want to create your own dataset please see the example GPSFM\Preprocessing\ProcessScript.m script that generates our format (Example pro.mat) from the standard format (Example.mat). 
projectivePipeline.m- script to run the pipeline. All the datasets from the paper are availble in the folder "DataSet Proj" and can be run using this script.
projectivePipelineSelfCalibrate.m- pipeline that performs also self calibration
runAll.m - script to run all the datasets in the paper. The script generates the results table and saves it to "table.csv" and also saves all the reconstructions after self calibration to the folder GPSFM\reconstructions



```
Important variables:

```
pointMatchesInliers- a matrix that contains the number of inliers between each pair of frames
FN- the multi-view fundamental matrix (measurements matrix)
M- tracks matrix 

```

## Acknowledgement 
We use the following 3rdparties code:
3rdparty\fromPPSFM - evaluation and self calibration code from "Practical projective structure
from motion (p2sfm)" (ICCV 2017). Their self calibration code implements the paper "Autocalibration via rank-constrained estimation of the
absolute quadric" (CVPR 2007), and uses the libraries "GloptiPoly 3" and "SeDuMi 1.3"

3rdparty\vgg_code - imlementation of basic functions from the book "Multiple View Geometry in Computer Vision" (2004) downloaded from:
https://www.robots.ox.ac.uk/~vgg/hzbook/code/

## Contact 
For any query, contact : 
Yoni Kasten, Amnon Geifman 
Weizmann Institute of Science
{yoni.kasten,amnon.geifman}@weizmann.ac.il

## License
   This software is provided under the provisions of the Lesser GNU Public License (LGPL). 
   see: http://www.gnu.org/copyleft/lesser.html.

   This software can be used only for research purposes, you should cite
   the aforementioned papers in any resulting publication.

   The Software is provided "as is", without warranty of any kind.




## Version History


* Version 1.0 (2019-06-06)
   Initial Release
