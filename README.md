# Polarimetry
Matlab codes for bundle calibration of a Mueller imaging polarimeter, and for Mueller polarimetric imaging, as presented in [1].

## Demos

The following demo files are provided: 

- `main_calibration.m` : demo of the calibration of a Mueller imaging polarimeter. Calibrates the zeros of the PSG polarizer, PSA and PSG retarders, plus the delays of the retarders.

- `main_polarimetry.m` : demo of Mueller imaging polarimetry, using the estimated calibration parameters.


## Images format

For the calibration codes to work, the following folder architecture should be used (see example in the Reflexion folder):

- the calibration data folder should contain N subfolders entitled with the wavelength in nanometers.
- each subfolder should contain (number of polarizer angles) x (number of PSG azimuts) x (number of PSA azimuts) images. The name of each image : angle1_angle2_angle3.format where angle1 is the PSA polarizer angle, angle2 is the PSG retarder azimut, and angle3 is the PSA retarder azimut. Angles in degrees should be rounded to the closest integer.

For the polarimetry codes to work, the following folder architecture should be used (see example in the Smileys folder):

- the calibration data folder should contain at least one subfolder entitled with the wavelength in nanometers.
- each subfolder should contain (number of PSG azimuts) x (number of PSA azimuts) images. The name of each image : angle2_angle3.format where angle2 is the PSG retarder azimut, and angle3 is the PSA retarder azimut. Angles in degrees should be rounded to the closest integer.


## References

[1] Yvain Quéau, Florian Leporcq, Ayman Alfalou. Design and simplified calibration of a Mueller imaging polarimeter for material classification. Opt.Lett., 2018, 43 (20), pp.4941-4944. DOI: 10.1364/OL.43.004941ff. Full text: https://hal.archives-ouvertes.fr/hal-01937915

Author of codes: Yvain Quéau, GREYC CNRS, yvain.queau@ensicaen.fr



