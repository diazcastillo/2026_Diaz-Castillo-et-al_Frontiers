This repository contains the R code necessary to replicate all analyses and figures presented in the publication: 
Diaz-Castillo C, Aguiar SR, Chamorro-Garcia R. Preconception exposures of female mice to a panel of metabolic disruptors induce sexually dimorphic metabolic perturbations in their offspring. Front Endocrinol. 2026;17. doi:10.3389/fendo.2026.1787973 (https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2026.1787973/full)

Prior to executing this code, the user must download all Supplementary files associated with this publication into a directory of their choice. 

The code will subsequently create three subdirectories and generate a single RData file.
  - The “temp” subdirectory will be utilized to store files necessary for advancing the analyses and will be deleted during the final step of the code.
  - The “File confirmations” subdirectory contains two additional subdirectories: “Supplementary files” and “Figures” The “Supplementary files” directory houses copies of Supplementary files for this publication that contain analytical results, including Supplementary Data 21-27 and Supplementary Figures 1 and 2. The “Figures” directory houses copies of the main figures of the publication. These copies will enable the user to compare them with the corresponding files accompanying this publication.
  - The RData file “2026_Diaz_Castillo_et_al.RData” will be periodically saved to store all main objects resulting from the code’s progress. This will provide the user with greater control over inspecting the code’s progress and even modifying it if desired.

Executing this code will require several days or hours depending on the user setup.
Commented lines in specific sections will inform the user of the steps that require the most time.
