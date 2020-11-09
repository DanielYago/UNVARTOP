# UNVARTOP
> Matlab code for 2D Topology Optimization using UNVARTOP method (for educational purposes)

<!-- BADGES -->
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/DanielYago/UNVARTOP/graphs/commit-activity)
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/DanielYago/UNVARTOP.svg)](https://GitHub.com/DanielYago/UNVARTOP/releases/)
[![GitHub issues](https://img.shields.io/github/issues/DanielYago/UNVARTOP.svg)](https://GitHub.com/DanielYago/UNVARTOP/issues/)
[![GitHub issues-closed](https://img.shields.io/github/issues-closed/DanielYago/UNVARTOP.svg)](https://GitHub.com/DanielYago/UNVARTOP/issues?q=is%3Aissue+is%3Aclosed)
[![GitHub contributors](https://img.shields.io/github/contributors/DanielYago/UNVARTOP.svg)](https://GitHub.com/DanielYago/UNVARTOP/graphs/contributors/)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html) 
[![DOI:10.1007/s00158-020-02722-0](https://zenodo.org/badge/DOI/10.1007/s00158-020-02722-0.svg)](https://doi.org/10.1007/s00158-020-02722-0)<!--
[![Github all releases](https://img.shields.io/github/downloads/DanielYago/UNVARTOP/total.svg)](https://GitHub.com/DanielYago/UNVARTOP/releases/)>
[![GitHub forks](https://img.shields.io/github/forks/DanielYago/UNVARTOP.svg?style=social&label=Fork&maxAge=2592000)](https://GitHub.com/DanielYago/UNVARTOP/network/)
[![GitHub watchers](https://img.shields.io/github/watchers/DanielYago/UNVARTOP.svg?style=social&label=Watch&maxAge=2592000)](https://GitHub.com/DanielYago/UNVARTOP/watchers/)
-->

<!-- TABLE OF CONTENTS -->

## Table of Contents

- [About the project](#about-the-project)
- [Getting Started](#getting-started)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
- [Executing program](#executing-program)
- [Features](#features)
- [Extensions](#extensions)
- [Contributing](#contributing)
- [Version History](#version-history)
- [Authors](#authors)
- [Support](#support)
- [License](#license)
- [Project Status](#project-status)
- [Sources](#sources)

<!-- ABOUT THE PROJECT -->
## About the project

This is the code for the paper: D. Yago, J. Cante, O. Lloberas-Valls, and J. Oliver. Topology Optimization using the UNsmooth VARiational Topology OPtimization (UNVARTOP) method: an educational implementation in Matlab. *Structural and Multidisciplinary Optimization*, 2020. [[10.1007/s00158-020-02722-0](https://doi.org/10.1007/s00158-020-02722-0)] 

The *unsmooth variational topology optimization* (UNVARTOP) approach is substantiated in a *bi-material setting*, via the *ersatz material* approach, defined by means of the *characteristic function*. Furthermore, sharp boundaries are obtained thanks to the use of a *discrimination function*, from which the *characteristic function* is computed. 

The optimal topology corresponds to the solution of a *closed-form algebraic system*, obtained from the topology optimization problem and involving the *relaxed topological derivative* (RTD) of the cost function. The resultant sensitivity of the cost function is regularized via a *Laplacian smoothing*, providing mesh-size control. The reference pseudo-time, imposed with the constraint equation and a *bisection algorithm*, is iteratively increased throughout the optimization, obtaining then intermediate converged optimal topologies, the so-called *incremental time-advancing scheme*. This approach provides not only the last optimal solution, but a set of optimal topologies for different volume percentages in a low number of iterations.

<!-- GETTING STARTED -->
## Getting Started

### Dependencies

* You must have Matlab (2017a or later) installed in your OS.
* Tested in Matlab 2018b (Windows 10) and 2019a (Debian 9).

### Installation

1. Clone this repository to your local machine
```shell
$ git clone https://github.com/DanielYago/UNVARTOP.git
```
2. Add directory to Matlab's path
```matlab
> install_UNVARTOP
```

<!-- EXECUTING PROGRAM -->
## Executing program

* Let us now run the following topology optimization:

  <p align="center">
	<img src="/Images/Cantilever_domain.png" alt="Domain" width="400">
	<!--![Domain](/Images/Cantilever_domain.png)-->
  </p>
  
  **Figure 1:** Cantilever beam: topology optimization domain and boundary conditions
  
  where a vertical load is applied at the bottom-right corner of the domain and the displacement is prescribed on the left side of it.
  
* According to the article, this numerical example can be carried out by evoking
  ```matlab
  UNVARTOP_2D_compliance (100,50,10,0,0.5,0,0.5)
  ```
  for `nelx` and `nely` equal to 100 and 50, respectively. The optimization will be performed in 10 equally-spaced time-steps, from   `t=0` to `t=0.5`. The dimensionless regularization parameter is set to 0.5 . 

* From the algorithm, the following graphs must be obtained

  <p align="center">
	<img src="/Images/Cantilever_topology.png" alt="Topology" width="400"> <img src="/Images/Cantilever_cost.png" alt="Cost_function" width="400">
	<!--![Topology](/Images/Cantilever_topology.png) ![Cost_function](/Images/Cantilever_cost.png)-->
  </p>
  **Figure 2 and 3:** Cantilever beam: optimal topology for `t=0.5` and evolution of the cost function, respectively.
  
* In addition, a GUI is generated with the topology evolution

  <p align="center">
	<img src="/Images/Cantilever_GUI.png" alt="GUI" width="450">
	<!--![GUI](/Images/Cantilever_GUI.png)-->
  </p>
  

 **Figure 4:** Cantilever beam: graphical user interface.

  where the deformed topology can be animated throughout the optimization by means of the buttons. For further information, please read the paper.

* Finally, the next animation can be generated from the set of time-steps

<p align="center">
   <img src="/Images/Cantilever_gif.gif" alt="Gif" width="600">
  <!--![Gif](/Images/Cantilever_gif.gif)-->
</p>

 **Figure 5:** Cantilever beam: topology optimization animation.

<!-- FEATURES -->
## Features

The algorithm is adapted in order to optimize:
 - [Minimum mean compliance](/Source_codes/Topology_optimization_problems/UNVARTOP_2D_compliance.m)
 - [Multi-load mean compliance](/Source_codes/Topology_optimization_problems/UNVARTOP_2D_multiload.m)
 - [Compliant mechanisms](/Source_codes/Topology_optimization_problems/UNVARTOP_2D_complmechanism.m)

Additionally, [UNVARTOP](/Source_codes/UNVARTOP_2D.m) combines the three topology optimization problems in a single code.

<!-- EXTENSIONS -->
## Extensions

Some of the extensions detailed in the paper are also implemented:
 - [Regula-falsi method](/Source_codes/Extensions/UNVARTOP_2D_compliance_RF.m)
 - [Anderson-Björck with Illinois algorithm](/Source_codes/Extensions/UNVARTOP_2D_compliance_AB.m)
 - [Augmented Lagrangian method](/Source_codes/Extensions/UNVARTOP_AL_2D.m)

Other numerical examples can be found in [the article](https://doi.org/10.1007/s00158-020-02722-0) or in [Examples](/Source_codes/Examples/). 

<!-- CONTRIBUTING -->
## Contributing

1. Fork or clone the repository: `$ git clone https://github.com/DanielYago/UNVARTOP.git`.
2. Create a new branch: `$ git checkout -b name_for_new_branch`.
3. Make changes and test them. Add new files to staging area: `$ git add name_of_new_file`.
4. Commit the changes: `$ git commit -am 'Add some comments'`.
5. Push changes to the new branch: `$ git push origin name_for_new_branch`.
6. Submit Pull Request with comprehensive description of changes

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on the process for submitting pull requests to us.

<!-- VERSION HISTORY -->
## Version History

* 1.1
	* Reviewed Release
    <!--* Various bug fixes and optimizations --->
* 1.0 
    * Submission Release
	* See [commit change](https://github.com/DanielYago/UNVARTOP/commits/master) or See [release history](https://github.com/DanielYago/UNVARTOP/releases)

<!-- AUTHORS -->
## Authors

* [Daniel Yago](https://github.com/DanielYago)
* [Juan Cante](https://github.com/jcante)
* [Oriol LLoberas-Valls](https://github.com/olloberas)
* Javier Oliver

<!-- SUPPORT -->
## Support

In case you find any errors or you have any suggestion, please send an email to `daniel.yago@upc.edu` or fork the repository and create a pull request.

<!-- LICENSE -->

## License

This project is licensed under the **[gpl-3.0](https://www.gnu.org/licenses/gpl-3.0.html)** License - see the [LICENSE](LICENSE) file for more details.

<!-- PROJECT STATUS -->
## Project Status

03/2020 - This project is part of the first author's PhD thesis, which is still in development in the field of topology optimization.

<!-- SOURCES -->
## Sources

[1] J. Oliver, D. Yago, J. Cante, and O. Lloberas-Valls. Variational approach to relaxed topological optimization: Closed form solutions for structural problems in a sequential pseudo-time framework. *Computer Methods in Applied Mechanics and Engineering*, 355:779-819, oct 2019. doi: [10.1016/j.cma.2019.06.038](https://doi.org/10.1016/j.cma.2019.06.038).
[2] D. Yago, J. Cante, O. Lloberas-Valls, and J. Oliver. Topology optimization of thermal problems in a nonsmooth variational setting: closed-form optimality criteria. *Computational Mechanics*, 66:259-286, june 2020. doi: [10.1007/s00466-020-01850-0](https://doi.org/10.1007/s00466-020-01850-0). 
[3] D. Yago, J. Cante, O. Lloberas-Valls, and J. Oliver. Topology Optimization using the UNsmooth VARiational Topology OPtimization (UNVARTOP) method: an educational implementation in Matlab. *Structural and Multidisciplinary Optimization*, 2020. doi: [10.1007/s00158-020-02722-0](https://doi.org/10.1007/s00158-020-02722-0). 