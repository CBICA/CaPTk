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
- This code (excluding dependent libraries) is governed by the license provided in http://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.
- The minimum recommended resolution is 1200x1024. We have seen some visualization issues with high DPI (>2K) screens and bug reports related to it will be appreciated.

## Download Latest Release (1.7.0)

By downloading CaPTk, you agree to our [License](./LICENSE).

| Platform (x64) | Link                                             |
|:--------------:|:------------------------------------------------:|
| Windows        | https://www.nitrc.org/frs/downloadlink.php/11203 |
| Linux          | https://www.nitrc.org/frs/downloadlink.php/11226 |
| macOS          | https://www.nitrc.org/frs/downloadlink.php/11212 |
| Archive        | https://www.nitrc.org/frs/?group_id=1059         |

## Frequently Asked Questions (FAQ)

### How do I run CaPTk Applications from the Command Line?

- The full list of command line applications available is shown by running the following command (pretty much every application is available via the command line):
```bash
${CaPTk_InstallDir}/bin/CaPTk -h
```

- An exemplary scenario to run Applications from the command line (the example shown is to get verbose help):

| Platform (x64) |   Functionality   |                              Successful Installation                             | Build From Source (after invoking "make install") | FUSE Issues (after invoking "--appimage-extract") |
|:--------------:|:-----------------:|:--------------------------------------------------------------------------------:|:-------------------------------------------------:|:-------------------------------------------------:|
|     Windows    |      Generic      |            ${CaPTk_InstallDir}/bin/${ApplicationName}.exe -h           | ${CaPTk_InstallDir}/bin/${ApplicationName}.exe -h |                        N.A.                       |
|      Linux     |      Generic      |              ${CaPTk_InstallDir}/bin/captk ${CaPTk_InstallDir}/bin/${ApplicationName}.cwl -h             |   ${CaPTk_InstallDir}/bin/${ApplicationName} -h   |      ~/CaPTk/${version}/${ApplicationName} -h     |
|      macOS     |      Generic      | ~/Applications/CaPTk_${version}.app/Contents/Resources/bin/${ApplicationName} -h |   ${CaPTk_InstallDir}/bin/${ApplicationName} -h   |                        N.A.                       |
|     Windows    | FeatureExtraction |            ${CaPTk_InstallDir}/bin/FeatureExtraction.exe -h            |  ${CaPTk_InstallDir}/bin/FeatureExtraction.exe -h |                        N.A.                       |
|      Linux     | FeatureExtraction |              ${CaPTk_InstallDir}/bin/captk FeatureExtraction.cwl -h              |    ${CaPTk_InstallDir}/bin/FeatureExtraction -h   |      ~/CaPTk/${version}/FeatureExtraction -h      |
|      macOS     | FeatureExtraction |  ~/Applications/CaPTk_${version}.app/Contents/Resources/bin/FeatureExtraction -h |    ${CaPTk_InstallDir}/bin/FeatureExtraction -h   |                        N.A.                       |

- There are detailed examples of individual command line usage in the [How To](https://cbica.github.io/CaPTk/How_To_Guides.html) section of the documentation.

### What are the OpenGL requirements to run CaPTk?

- If CaPTk is unable to load images or you receive the error about minimum OpenGL version wasn't found, please update your display drivers in order to have **OpenGL version 3.2 or above**. Some useful resources:
  - OpenGL update for Ubuntu [[ref](https://www.phoronix.com/scan.php?page=news_item&px=Ubuntu-16.04-OI-Intel-GL-4.2)]: `sudo apt-add-repository ppa:oibaf/graphics-drivers && sudo apt-get update && sudo apt-get dist-upgrade`
  - https://community.khronos.org/t/how-to-update-opengl/75314
  - https://ubuntuforums.org/showthread.php?t=2326268
  - https://www.techwalla.com/articles/how-to-update-opengl-drivers

### What if I am having issues with the Linux Installer?

- If the installer successfully finishes and you are not able to run CaPTk due to FUSE issues, please extract the installer using the following command to extract the contents of the AppImage onto the hard drive: 
```bash
user@pc:~# ~/CaPTk/${version}/captk --appimage-extract
```
  - This will extract the package to the path `~/CaPTk/${version}/squashfs-root/usr/` with the binaries present in `~/CaPTk/${version}/squashfs-root/usr/bin`.
  - To run **CaPTk** from the command line in this manner, the user would need to make the following additions:
    - Add the path `~/CaPTk/${version}/squashfs-root/usr/lib` to their environment variable **PATH**: `export PATH=~/CaPTk/${version}/squashfs-root/usr/bin:$PATH`
    - Add the path `~/CaPTk/${version}/squashfs-root/usr/lib` to their environment variable **LD_LIBRARY_PATH**: `export LD_LIBRARY_PATH=~/CaPTk/${version}/squashfs-root/usr/lib:$LD_LIBRARY_PATH`
    - This will need to be present in the user's shell start up file `.bashrc`, `.cshrc`, `.profile`, `.kshrc`, etc.` for the settings to be persistent across login sessions.
    - All the command line applications should now be available for use. For example: 
    ```bash
    FeatureExtraction \ # the executable should already be in the $PATH
      -n AAC0_timestamp \ # this is the subject ID for this run 
      -i /usr/path/T1.nii.gz,/usr/path/T2.nii.gz -t T1,T2 \ # the co-registered input image(s) and their respective modalities
      -m /user/path/mask.nii.gz \ # the mask file (co-registered with the input image(s)
      -r 2,4,5 -l ED,EN,NE \ # the labels from the mask file from where features need to be extracted and their respective label identifier(s)
      -p /usr/path/features.csv \ # the parameter file to use 
      -o /usr/path/output.csv # the path to write the output 
    ```
- Currently, we support all distributions newer than Ubuntu 16.04.

### What is the compatibility of various CaPTk installers?
  
| Platform (x64) |                         Tested                         |       Unsupported      |
|:--------------:|:------------------------------------------------------:|:----------------------:|
|     Windows    |                        7, 8, 10                        |        XP, Vista       |
|      Linux     | Ubuntu 16.04, 18.04; Debian 9, CentOS 7 (source build) | Ubuntu 14.04; CentOS 6 |
|      macOS     |                      10.13, 10.14                      |          10.12         |

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
