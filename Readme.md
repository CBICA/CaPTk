# CaPTk:  Cancer Imaging Phenomics Toolkit 

<p align="center">
    <img src="https://raw.githubusercontent.com/CBICA/CaPTk/master/docs_sources/images/0_tutorial_1_mainUI_resized.png" />
    <br></br>
    <a href="https://dev.azure.com/CBICA/CaPTk/_build?definitionId=2" alt="Build Status"><img src="https://dev.azure.com/CBICA/CaPTk/_apis/build/status/CBICA.CaPTk?branchName=master" /></a>
    <a href="https://hub.docker.com/r/cbica/captk/builds" alt="Build Status"><img src="https://img.shields.io/docker/cloud/automated/cbica/captk"></a>
    <a href="https://github.com/CBICA/CaPTk/issues" alt="Issues"><img src="https://img.shields.io/github/issues/CBICA/CaPTk.svg" /></a>
    <a href="https://github.com/CBICA/CaPTk/issues" alt="Issues"><img src="https://img.shields.io/github/issues-closed/CBICA/CaPTk.svg" /></a>
    <a href="https://doi.org/10.1117/1.JMI.5.1.011018" alt="Citation"><img src="https://img.shields.io/badge/cite-citation-blue" /></a>
    <img src="https://img.shields.io/badge/language-c%2B%2B11-blue.svg" />
</p>

CaPTk is a software platform, written in C++, for analysis of radiographic images of cancer, currently focusing on brain, breast, and lung cancer. CaPTk integrates advanced, validated tools performing various aspects of medical image analysis, that have been developed in the context of active clinical research studies and collaborations toward addressing real clinical needs. With emphasis given in its use as a very lightweight and efficient viewer, and with no prerequisites for substantial computational background, CaPTk aims to facilitate the swift translation of advanced computational algorithms into routine clinical quantification, analysis, decision making, and reporting workflow.

Its long-term goal is to provide widely used technology that leverages the value of advanced imaging analytics in cancer prediction, diagnosis, and prognosis, as well as in better understanding the biological mechanisms of cancer development.

CaPTk is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania.

For a full list of applications and functionalities, please see https://cbica.github.io/CaPTk/How_To_Guides.html

For more details, please visit us at https://www.cbica.upenn.edu/captk

For project documentation and how-to guides, please visit https://cbica.github.io/CaPTk/

For issues, please visit https://github.com/cbica/captk/issues

## Supporting Grant
This work is in part supported by the grant U24-CA189523, awarded by the National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research (NIH/NCI/ITCR).

## Disclaimer
- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- This code (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.

## Downloads

By downloading CaPTk, you agree to our [License](./LICENSE). You can review Installation Instructions [here](https://cbica.github.io/CaPTk/Installation.html).

## Latest Stable (1.7.6)

| Platform (x64) | Link                                             |
|:--------------:|:------------------------------------------------:|
| Windows        | https://www.nitrc.org/frs/downloadlink.php/11652 |
| Linux          | https://www.nitrc.org/frs/downloadlink.php/11651 |
| macOS          | https://www.nitrc.org/frs/downloadlink.php/11650 |
| Archive        | https://www.nitrc.org/frs/?group_id=1059         |

## Test Build (1.8.0.Alpha2)

| Platform (x64) | Link                                             |
|:--------------:|:------------------------------------------------:|
| Windows        | https://www.nitrc.org/frs/downloadlink.php/11777 |
| Linux          | https://www.nitrc.org/frs/downloadlink.php/11778 |
| Linux (CentOS7)| https://www.nitrc.org/frs/downloadlink.php/11779 |
| macOS          | https://www.nitrc.org/frs/downloadlink.php/11780 |

## Development Builds

These are UNTESTED development builds from the latest master. Use at your own risk.

| Platform (x64) | Link                                             |
|:--------------:|:------------------------------------------------:|
| Windows        | https://www.nitrc.org/frs/downloadlink.php/11516 |
| Linux          | https://www.nitrc.org/frs/downloadlink.php/11517 |
| Linux (CentOS7)| https://www.nitrc.org/frs/downloadlink.php/11613 |
| macOS          | https://www.nitrc.org/frs/downloadlink.php/11518 |

## [Frequently Asked Questions (FAQ)](https://cbica.github.io/CaPTk/Getting_Started.html#gs_FAQ)

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
