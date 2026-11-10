# Sources and observation provenance

The implementation uses numerical controls, not an observational calibration. These sources set the context and identify effects that a successful full model must include.

- [Torii et al. (2024), Global N-body Simulation of Gap Edge Structures](https://arxiv.org/html/2403.03012v1): resolved collective dynamics and a collision-driven vertical-wall mechanism. The study uses large N-body calculations. Its discussion gives variable edge distances rather than one fixed impact parameter. This motivates the radial sweep; it does not validate our smaller boxes.
- [Torii et al. (2025), effect of eccentricity and inclination](https://arxiv.org/abs/2504.10989): eccentric forcing can localize strong collisions and vertical structures. The paper reports a Hill-scaled wall height near 0.2 Hill radii. It favors collision-driven walls over direct vertical forcing. Our separate experiments do not independently reproduce that full result.
- [Seiß et al. (2010), The vertical structure of the Daphnis wakes](https://meetingorganizer.copernicus.org/EPSC2010/EPSC2010-701.pdf): an earlier study of collisional vertical excitation and satellite inclination.
- [Tiscareno et al. (2008), Propeller Population in Saturn's A Ring, section 6](https://mmhedman.github.io/papers_published/MST_proppop_aj_2008.pdf): discusses a brightness increase without an increase in mass density and the possible disruption of opaque self-gravity wakes. This is a candidate explanation, not a measurement of the images in this task.
- [Rein and Liu (2012), REBOUND](https://www.aanda.org/articles/aa/pdf/2012/01/aa18085-11.pdf), and the [official shearing-sheet example](https://rebound.hanno-rein.de/c_examples/shearing_sheet/): the numerical engine and the existing ring-box method.
- [Mondino-Llermanos and Salo (2022)](https://academic.oup.com/mnras/article/513/4/4711/6575940): explains why dense or sustained contacts can require soft-sphere treatment. Binary hard-sphere tests do not validate that regime.
- [JPL satellite physical parameters](https://ssd.jpl.nasa.gov/sats/phys_par/): reference for satellite properties. Model constants remain explicit nominal inputs in `encounter.py`.

## The two supplied images

[PIA08319, as supplied on Commons](https://commons.wikimedia.org/wiki/File:PIA08319_Daphnis_in_Keeler_Gap.jpg) was acquired on October 27, 2006. The caption gives a lit-side view about 25° below the ring plane, phase angle 36°, and original scale 2 km/pixel. The uploaded image was enlarged by a factor of two and sharpened. It is not a calibrated radiance product; enlarged pixels do not provide new resolution. Credit: NASA/JPL/Space Science Institute.

[The supplied shadow image](https://commons.wikimedia.org/wiki/File:Daphnis_edge_wave_shadows.jpg) is **PIA11653**, not the different image titled PIA11579 “Edge Wave Shadows.” The supplied image was acquired on May 24, 2009, with a lit-side view about 57° below the plane, phase angle 121°, and a stated scale of 5 km/pixel. The caption reports structures 0.5–1.5 km high. It does not give the solar elevation needed to infer height from a measured shadow length. Credit: NASA/JPL/Space Science Institute.

[NASA's related Rippling Shadows caption](https://science.nasa.gov/resource/rippling-shadows/) provides the same reported height range. These heights concern visible structures, not the standard deviation of every particle's vertical position. A comparison must use the same height statistic, projection, optical threshold, and image resolution.
