# Scheid_RNS-IEEG-Biomarkers
Tools for biomarker analysis and figure generation used in the following publication: 

Scheid BH, Bernabei JM, Khambhati AN, Mouchtaris S, Jeschke J, Bassett DS, et al. Intracranial electroencephalographic biomarker predicts effective responsive neurostimulation for epilepsy prior to treatment. Epilepsia. 2022 Jan 7;epi.17163. 

### Requirements
- Matlab >= 2018a
- User account on [ieeg.org](http://ieeg.org)
- [ieeg cli](https://bitbucket.org/ieeg/ieeg/wiki/cli.md)

### Data Access
* Create an account on ieeg.org
* Download the [ieeg cli](https://bitbucket.org/ieeg/ieeg/wiki/cli.md)
* Add `ieeg` to your PATH variable as a link to ieeg.bin (usually added as /usr/local/bin/ieeg-cli-1.xx.xx/ieeg) 
* Use the following commands to access the data:

```
ieeg download-rec-objects --dataset UCSF_Biomarkers_Data UCSF_shared.mat
ieeg download-rec-objects --dataset NYU_Biomarkers_Data NYU_shared.mat
ieeg download-rec-objects --dataset Penn_Biomarkers_Data Penn_shared.mat
```

### Running the Pipeline
1. Add path to figure folder
2. Add path to data folder
3. Hit Run!

