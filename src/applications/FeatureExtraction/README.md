# Feature Extraction

## Command Line Usage

### Windows PowerShell / CMD:
```
FeatureExtraction.exe –o C:\path\to\output\csv –b C:\path\to\batch\csv –p C:\path\to\parameters\csv [– d 1] [– vc 1]
```

### Linux/macOS Terminal:
```
FeatureExtraction –o /path/to/output/csv –b /path/to/batch/csv –p /path/to/parameters/csv [– d 1] [– vc 1]
```

## Locating the executable

### Windows

Under the CaPTk installation tree, you will see the folder */bin/*; the FeatureExtraction executable will be there.

### macOS

Under the *Applications* folder, ```cd``` into ```CaPTk.app/Contents/Resources/bin``` and the FeatureExtraction will be there.

### Linux

The Linux package is an [AppImage](https://appimage.org/) and the FeatureExtraction executable can be access using the following [cwl](https://www.commonwl.org/)-enabled command (assuming you have the CaPTk AppImage in your PATH):

```captk featureextraction $your_command_goes_here```