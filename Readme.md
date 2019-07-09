# CaPTk:  Cancer Imaging Phenomics Toolkit 

<p align="center">
    <img src="https://www.med.upenn.edu/cbica/assets/user-content/images/captk/baseScreenshot.png" />
    <br></br>
    <a href="https://dev.azure.com/CBICA/CaPTk/_build?definitionId=2" alt="Build Status"><img src="https://dev.azure.com/CBICA/CaPTk/_apis/build/status/CBICA.CaPTk?branchName=master" /></a>
    <a href="https://github.com/CBICA/CaPTk/issues" alt="Issues"><img src="https://img.shields.io/github/issues/CBICA/CaPTk.svg" /></a>
    <a href="https://github.com/CBICA/CaPTk/issues" alt="Issues"><img src="https://img.shields.io/github/issues-closed/CBICA/CaPTk.svg" /></a>
    <img src="https://img.shields.io/badge/language-c%2B%2B11-blue.svg" />
</p>

CaPTk is a software platform, written in C++, for analysis of radiographic images of cancer, currently focusing on brain, breast, and lung cancer. CaPTk integrates advanced, validated tools performing various aspects of medical image analysis, that have been developed in the context of active clinical research studies and collaborations toward addressing real clinical needs. With emphasis given in its use as a very lightweight and efficient viewer, and with no prerequisites for substantial computational background, CaPTk aims to facilitate the swift translation of advanced computational algorithms into routine clinical quantification, analysis, decision making, and reporting workflow.

Its long-term goal is to provide widely used technology that leverages the value of advanced imaging analytics in cancer prediction, diagnosis, and prognosis, as well as in better understanding the biological mechanisms of cancer development.

CaPTk is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania.

For more details, please visit us at https://www.cbica.upenn.edu/captk

For project documentation and how-to guides, please visit https://cbica.github.io/CaPTk/

For issues, please visit https://github.com/cbica/captk/issues

## Supporting Grant
This work is in part supported by the grant U24-CA189523, awarded by the National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research (NIH/NCI/ITCR).

## Disclaimer
- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- This code (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.

## Download Latest Release (1.7.1)

By downloading CaPTk, you agree to our [License](./LICENSE).

| Platform (x64) | Link                                             |
|:--------------:|:------------------------------------------------:|
| Windows        | https://www.nitrc.org/frs/downloadlink.php/11479 |
| Linux          | https://www.nitrc.org/frs/downloadlink.php/11480 |
| macOS          | https://www.nitrc.org/frs/downloadlink.php/11481 |
| Archive        | https://www.nitrc.org/frs/?group_id=1059         |

## [Frequently Asked Questions (FAQ)](https://sarthakpati.github.io/CaPTk/gs_FAQ.html)

## Contact
For more information, please contact <a href="mailto:software@cbica.upenn.edu">CBICA Software</a>.

## GitHub Distribution

We currently provide only our tagged versions of the code via GitHub. Check the "tags" using your favorite Git client after cloning our repository. The analogous commands are as follows:

```bash
git clone https://github.com/cbica/captk.git
latesttag=$(git describe --tags)
echo checking out ${latesttag}
git checkout ${latesttag}
```
