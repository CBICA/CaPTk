# Frequently Asked Questions (FAQ)

<details>
  <summary>Which platforms are supported by CaPTk installers?</summary>
  
| Platform (x64) |                         Tested                         |       Unsupported      |
|:--------------:|:------------------------------------------------------:|:----------------------:|
|     Windows    |                        7, 8, 10                        |        XP, Vista       |
|      Linux     | Ubuntu 16.04, 18.04; Debian 9, CentOS 7 (source build) | Ubuntu 14.04; CentOS 6 |
|      macOS     |                      10.13, 10.14                      |          10.12         |

</details>

<details>
  <summary>How do I run CaPTk Applications from the Command Line?</summary>

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
</details>

<details>
  <summary>What are the OpenGL requirements to run CaPTk?</summary>

- If CaPTk is unable to load images or you receive the error about minimum OpenGL version wasn't found, please update your display drivers in order to have **OpenGL version 3.2 or above**. Some useful resources:
  - OpenGL update for Ubuntu [[ref](https://www.phoronix.com/scan.php?page=news_item&px=Ubuntu-16.04-OI-Intel-GL-4.2)]:
  ```bash
  sudo apt-add-repository ppa:oibaf/graphics-drivers && sudo apt-get update && sudo apt-get dist-upgrade
  ```
  - https://community.khronos.org/t/how-to-update-opengl/75314
  - https://ubuntuforums.org/showthread.php?t=2326268
  - https://www.techwalla.com/articles/how-to-update-opengl-drivers
</details>

<details>
  <summary>What if I am having issues with the Linux Installer?</summary>

### FUSE Issues

If the installer successfully finishes and you are not able to run CaPTk due to FUSE issues, please extract the installer using the following command to extract the contents of the AppImage onto the hard drive: 

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
    
### GLIBCXX or GLIBC Issues

- This is happening because the Linux binaries are compiled using GCC 4.8.5 on Ubuntu 16.04; therefore you will need to update GCC 4.8.5 or newer in order to get CaPTk to work. Some examples are shown:

  - Ubuntu [[ref](https://askubuntu.com/a/581497)]:

  ```bash
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test
  sudo apt-get update
  sudo apt-get install gcc-4.9 g++-4.9
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9
  ```  

  - CentOS [[ref](https://www.softwarecollections.org/en/scls/user/rhscl/?search=devtoolset&policy=&repo=&order_by=-create_date&per_page=10)]:

  ```bash
  sudo yum install centos-release-scl
  sudo yum install devtoolset-6
  ```  

- You can still build CaPTk from source using the instructions in [Technical Reference](https://cbica.github.io/CaPTk/Technical_Reference.html).

</details>

<details>
  <summary>What if I am having an issue not listed here?</summary>
  
Please [open a new issue](https://github.com/CBICA/CaPTk/issues/new?assignees=&labels=&template=bug-report.md&title=) with us and we will do our best to resolve it.

</details>
