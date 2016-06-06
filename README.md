# VBMplots
Shiny app for visualization voxel-based analysis results.
## How to Use
To use this shiny app you should have results of your analysis as nifti file in MNI space (197x233x189), where per every voxel your stored result value in corresponding coordinate.

In **data** folder create new directory **results** and per every nifti image of your results make separate folder in put image there:

```
-data\
--results\
---APOE4\
----apoe4_p_value.nii.gz \
```

By default VBMplots app asume that image contain p-values therefore apply -log10 transformation. To skip transformation step add in your **server.R** file to "Read_region_data" and "Read_data" function additional parameter **data_type='not p-value'**   

## How to Run

1) install shiny library : `install.packages('shiny')`

2) import library: `library(shiny)`

3) Navigate to parent directory of **VBMplots** folder and run app: `runApp('VBMplots')`



## Example
Check the example: [AD SNPs plots](http://www.roshchupkin.com/adsnps/)

## Advanced user

You can easily change the atlas for brain region labeling for you own one (in this case the result image could be in any space and have any resolution).
Use **Atlas.txt** table as example to make and replace by your own table linked to your Atlas. Then replace  **Atlas.nii.gz** image with your atlas, where voxel value coded structure from your Atlas table.