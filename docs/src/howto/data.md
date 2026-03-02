# Obtaining input data

The minimal input data required to run the model will have been downloaded automatically
from Zenodo during installation. If you require more data, such as additional stellar
spectra or opacities, these can also be obtained using the `get_data` script in the AGNI
root directory. To see how to use this script, run it without arguments:
```bash
./src/get_data.sh
```

## Spectral files

Opacities are contained within "spectral files". Use the table within
`res/spectral_files/reference.pdf` to decide which spectral files are best for you.

For example, if you wanted to get the spectral file "Honeyside48" you would run:
```bash
./src/get_data.sh anyspec Honeyside 48
```

Additional spectral files can also be downloaded directly from the
[PROTEUS community on Zenodo](https://zenodo.org/communities/proteus_framework/records?q&f=subject%3Aspectral_files&l=list&p=1&s=10&sort=newest).

!!! tip
    See [Spectral file does not exist](@ref) in the troubleshooting guide if a spectral
    file you downloaded cannot be found by AGNI.
