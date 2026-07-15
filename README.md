<!-- PROJECT SHIELDS -->
[![arXiv][arxiv-shield]][arxiv-url]
[![MIT License][license-shield]][license-url]
[![ReseachGate][researchgate-shield]][researchgate-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
[![Scholar][scholar-shield]][scholar-url]
<!-- [![finalpaper][finalpaper-shield]][finalpaper-url] -->
<!-- [![Webpage][webpage-shield]][webpage-url] -->

# A Generalized Plant Perspective on Linear-Convex Feedback Optimization

This repository contains the implementations from our paper

> F. Jakob and A. Iannelli. "A Generalized Plant Perspective on Linear-Convex Feedback Optimization." arXiv preprint. arXiv:2606.14471 (2026). 

## Installation

This repository was developed and tested using MATLAB R2023b with the following dependencies:

- Simulink
- Robust Control Toolbox
- [IQClab](https://github.com/JoostVeenman/IQClab)

Standalone ODE solver routines will follow to replace the Simulink dependency.

### Installing IQClab

The correct commit of IQClab will be pulled as submodule upon cloning:

```
git clone  --recurse-submodules https://github.com/col-tasas/2026-iqc-feedback-opt.git
```

After cloning, go to the [IQClab installation file](https://github.com/JoostVeenman/IQClab/blob/main/IQClab_install.m) and correctly set the path of the toolbox.

IQClab can be used standalone through the MATLAB built-in SDP parser LMILab. However, several scripts rely on the alternative parser YALMIP in combination with the MOSEK solver. YALMIP can be downloaded [here](https://yalmip.github.io/download/), and an academic MOSEK license can be requested [here](https://www.mosek.com/products/academic-licenses/).


## Running the code

The `src` directory yields executable scripts to recreate the numeric experiments of the paper.

Note that the IQC synthesis is sensitive w.r.t. its numerical parameters. Careful tuning is recommended.


## Contact
🧑‍💻 Fabian Jakob

📧 [fabian.jakob@ist.uni-stuttgart.de](mailto:fabian.jakob@ist.uni-stuttgart.de)



[license-shield]: https://img.shields.io/badge/License-MIT-T?style=flat&color=blue
[license-url]: https://github.com/col-tasas/2026-iqc-feedback-opt/blob/main/LICENSE
[webpage-shield]: https://img.shields.io/badge/Webpage-Fabian%20Jakob-T?style=flat&logo=codementor&color=green
[webpage-url]: https://www.ist.uni-stuttgart.de/de/institut/team/Jakob-00004/
[arxiv-shield]: https://img.shields.io/badge/arXiv-2606.14471-t?style=flat&logo=arxiv&logoColor=white&color=red
[arxiv-url]: https://arxiv.org/abs/2606.14471
[finalpaper-shield]: https://img.shields.io/badge/IEEE-Paper-T?style=flat&color=red
[finalpaper-url]: https://google.com
[researchgate-shield]: https://img.shields.io/badge/ResearchGate-Fabian%20Jakob-T?style=flat&logo=researchgate&color=darkgreen
[researchgate-url]: https://www.researchgate.net/profile/Fabian-Jakob-4
[linkedin-shield]: https://img.shields.io/badge/Linkedin-Fabian%20Jakob-T?style=flat&logo=linkedin&logoColor=blue&color=blue
[linkedin-url]: https://de.linkedin.com/in/fabian-jakob
[scholar-shield]: https://img.shields.io/badge/Google-Scholar-T?style=flat&logo=googlescholar&color=blue
[scholar-url]: https://scholar.google.com/citations?user=WQsMJp0AAAAJ&hl=en
