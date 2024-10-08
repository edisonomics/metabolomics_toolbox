## Example workflows for metabolomics_toolbox
The same content, and additional material, are also available on NMRbox:
**`/public/metabolomics-toolbox/example_matlab_workflows`**

|Study Type|Dataset name|About the Study|This Dataset Contains|This Workflow Uses|
|-|-|-|-|-|
| `1D` `1H` `urine` `human`|[matlab_workflow1_complete_nan](https://github.com/edisonomics/metabolomics_toolbox/tree/master/examples/1D_serum/matlab_workflow1_complete_nan)|Urine NMR metabolomics study. This study was originally conducted by Olatomiwa Bifarin (J. Proteome Res., 2021)|[Spectra](https://github.com/edisonomics/metabolomics_toolbox/tree/master/examples/1D_serum/matlab_workflow1_complete_nan/data/spectra) [Workflow](https://github.com/edisonomics/metabolomics_toolbox/blob/master/examples/1D_serum/matlab_workflow1_complete_nan/scripts/matlab_workflow1_complete_nan.m)|`Load1D` `Setup1D` `displaypeak1D` `ref_spectra` `remove_region` `guide_align1D` `normcheck` `normalize` `varcheck` `scale` `nipalsPCA` `VisScores`|
||matlab_workflow1_summarized_nan*|Summarized version of `matlab_workflow1_complete_nan`. This version uses pre-processed data and requires less time to complete the workflow.|||

*: Only available on NMRbox due to the large file size
