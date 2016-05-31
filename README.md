# VBMplots
Shiny app for visualization voxel-based analysis results
## How to Use
To use this shiny app you should have results of your analysis as nifti file in MNI space (197x233x189), where per every voxel your stored result value in corespondent coordinate.

In **data** folder create new directory **results** and per every nifti image of your results make separate folder in put image there:

```
-data\
--results\
---APOE4\
----apoe4_p_value.nii.gz \
```

By default VBMplots app asume that image contain p-values therefore apply -log10 strasformation. To skip trasformation step add in your **server.R** file to "Read_region_data" and "Read_data" function additinal parameter **data_type='not p-value'**   

## Example
Check the example: [AD SNPs plots](http://www.roshchupkin.com/adsnps/)
