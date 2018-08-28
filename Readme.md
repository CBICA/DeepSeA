# DeepSeA
### ACKNOWLEDGEMENT
If you use this program in your research and/or publication(s), you must cite the paper below:
	[1] D. Wei, S. Weinstein, M.-K. Hsieh, L. Pantalone, M. Schnall, and D. Kontos, "Three-dimensional whole breast segmentation in sagittal MR images with dense depth field modeling and localized self-adaptation," in Medical Imaging: Image Processing, ser. Proc. SPIE, vol. 10133, 2017, pp. 1013314.

### FUNCTIONALITY
This package provides functionality for 3D whole breast segmentation in breast MRI. It supports both sagittal and axial breast MRI data. Currently it requires T1w nonfat-suppressed images for segmentation.

### ENVIRONMENTS
This package has been tested under:
1. Windows 7 Ultimate with Service Pack1, using MATLAB R2014b.
2. Windows 7 Enterprise with Service Pack1, using MATLAB R2014b.
3. Windows 10 Home Edition 17134.165, using MATLAB R2017a.

### USAGE
To start using the package, call
	
```matlab
MR_startup.m
```

which will add all modules into MATLAB paths so that user can call the pipeline and individual functions in the workspace.

### EXAMPLES
This package includes two examples with anonymized real subject data --- one for sagittal breast MRI and the other for axial breast MRI --- to illustrate basic usage.

Example #1:

```matlab 
segdata = wholeBreastSegment('demoData\AX_T1_BILAT\', 'Results\', 'LateralBounds', true)
```

This example loads axial breast MRI data (T1w nonfat-suppressed) from DICOM files in the folder ```.\demoData\AX_T1_BILAT\```, and writes the results to the folder ```.\Results\```. In addition, it sets the optional input 'LateralBounds' to true, dictating the program to define lateral bounding planes for the breasts on the outer sides of the body. The output ```segdata``` is a struct containing various segmentation results, e.g., 3D breast masks, breast size in cm^3, etc. Please refer to the help info of the function ```wholeBreastSegment``` for detailed description.

Example #2:

```matlab
segdata = wholeBreastSegment('demoData\SAG_T1_3D_BILAT\', 'Results\', 'BottomBound', true)
```
This example loads sagittal breast MRI data (T1w nonfat-suppressed) from DICOM files in the folder ```.\demoData\SAG_T1_3D_BILAT\```, and writes the results to the folder ```.\Results\```. In addition, it sets the optional input 'BottomBound' to true, dictating the program to define the inferior bounding planes for the breasts. The output ```segdata``` is a struct containing various segmentation results, e.g., 3D breast masks, breast size in cm^3, etc. Please refer to the help info of the function ```wholeBreastSegment``` for detailed description.
Example #2 takes much longer than Example #1 to finish due to the finer in-plane resolution of the sagittal demo data. However, you can also control the resolution used during processing by setting the ```VoxelSize``` optional input; please refer to the help info of the function ```wholeBreastSegment``` for detailed description.
	
### CONTACT
For bug regporting, etc., please email:
Dong.Wei@uphs.upenn.edu / Meng-Kang.Hsieh@uphs.upenn.edu
