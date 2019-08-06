# Binaries

Contains the binaries used by CI/the superbuild. These are automatically handled for you by cmake, no need to touch them.

## License Agreement

By downloading these binaries, you are assuming acceptance of the specific licenses being described in ../licenses/Combined.txt

## Qt Options

Extraction has been done using the relevant Qt installer from the web by enabling the following options for the chosen compiler:

For project documentation and how-to guides, please visit https://cbica.github.io/CaPTk/

For issues, please visit https://github.com/cbica/captk/issues

## Supporting Grant
This work is in part supported by the grant U24-CA189523, awarded by the National Institutes of Health / National Cancer Institute / Informatics Technology for Cancer Research (NIH/NCI/ITCR).

## Disclaimer
- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- This code (excluding dependent libraries) is governed by the license provided in https://www.med.upenn.edu/sbia/software-agreement.html unless otherwise specified.

## Downloads

### Latest Release (1.7.2)

By downloading CaPTk, you agree to our [License](./LICENSE).

| Platform (x64) | Link                                             |
|:--------------:|:------------------------------------------------:|
| Windows        | https://www.nitrc.org/frs/downloadlink.php/11503 |
| Linux          | https://www.nitrc.org/frs/downloadlink.php/11487 |
| macOS          | https://www.nitrc.org/frs/downloadlink.php/11486 |

## Older Releases

The entire archive of CaPTk's binaries in our NITRC Download page: https://www.nitrc.org/frs/?group_id=1059

## [Frequently Asked Questions (FAQ)](https://cbica.github.io/CaPTk/gs_FAQ.html)

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
