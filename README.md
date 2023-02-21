# Imaging-Glass-and-Mirrors
Code and data associated with the paper "Detection and Mapping of Specular Surfaces Using Multibounce Lidar Returns" by Henley et al. 2022"

## Code 
(see comments for more detailed descriptions):

### Raw Data Processing
  ***batchProcessMirrorMultiBounce.m***
  
  Driver script that batch processes raw photon count measurements (spot_XX.csv).  Outputs are data cubes as well as per-pixel maps of quanities like pixel-wise maps of detected pulse energy and peak timing bins for each laser scan direction (spot_XX.mat).
  
  ***processRawCounts.m***
  
  From raw photon count measurement file (spot_XX.csv), outputs data cube for which photon counts are binned by time-of-arrival and detector scan angle.
  
  ***datacubeStats.m***
  
  Using photon count data cube as input, outputs per-pixel maps of quantities such as per-pixel detected pulse energy and timing bin of detected peak.
    
### Spot Extraction
***processFrame.m***

Reads in frame information, detects visibile specularities and output relevant properties, including spot energy, time-of-flight, and angle of arrival.

***spotDetector.m***

Detects laser spots in data and extracts their angle of arrival and intensity (in photon counts).
 
***binDetector.m***

Extracts times-of-flight associated with all detected laser spots.

### Specular surface mapping with single-beam illumination

***mirrorGeometryScript.m***
Driver script to generate specular and diffuse point clouds from measurements.

***computeGeometryMirrorDisambiguation.m***
Compute 3D points from spot detections.  Includes disambiguation logic for scenes that contain transparent surfaces.

***specularPointCloudPlotter.m***
Plots the point cloud generated using mirrorGeometryScript.  Points are colored by reflection type.

### Specular surface mapping with multi-beam illumination
***processMultibeamData.m***
Driver script to implement specular surface mapping using multi-beam illumination.

***mirrorSourceFit.m***
Robustly localizaes the position of the mirrored source using a RANSAC algorithm.

### Miscellaneous
***naiveMirrorGeometryScript.m***
Driver script that generates point cloud from measurements assuming that all detections correspond to one-bounce light transport paths.       

***computeGeometryNaively.m***
Computes 3D points from detections.  Assumes all detected spots are one-bounce returns.

***naivePointCloudPlotter.m***
Plots the point cloud generated using naiveMirrorGeometryScript.m.
                       

               
                      
