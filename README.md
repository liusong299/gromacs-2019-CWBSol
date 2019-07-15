# Gromacs CWBSol
The original release of CWBSol method with Gromacs 2019

## Getting Started
The 3D reference interaction site model (CWBSol) is a powerful tool to study the thermodynamic and structural properties of liquids. However, for hydrophobic solutes, the inhomogeneity of the solvent density around them poses a great challenge to the 3DRISM-HI theory. To address this issue, we have previously introduced the hydrophobic-induced density inhomogeneity theory (HI) for purely hydrophobic solutes. To further consider the complex hydrophobic solutes containing partial charges, here we propose the D2MSA closure to incorporate the short-range and long-range interactions with the D2 closure and the mean spherical approximation, respectively. We demonstrate that our new theory can compute the solvent distributions around real hydrophobic solutes in water and complex organic solvents that agree well with the explicit solvent molecular dynamics simulations.

### Prerequisites

#### I. Package requirement
* [FFTW3](http://www.fftw.org) with double precision (required)
* [GROMACS 2019](http://www.gromacs.org) (required)

#### II.Installing

Introduction to building GROMACS with CWBSol method.These instructions pertain to building GROMACS 2019 with CWBSol method. You might also want to check the up-to-date installation instructions from original Gromacs installation.

----------------------------

1. Get the latest version of your C and C++ compilers.

2. Check that you have CMake version 3.4.3 or later.

3. Get and unpack the latest version of the GROMACS tarball.

4. Make a separate build directory and change to it.

5. Run "cmake" with the path to the source as an argument

6. Run "make", "make check", and "make install"

7. Source "GMXRC" to get access to GROMACS

Or, as a sequence of commands to execute:
```
   tar xfz gromacs-2019.tar.gz
   cd gromacs-2019
   mkdir build
   cd build
   cmake .. -DGMX_DOUBLE=ON -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
   make
   make check
   sudo make install
   source /usr/local/gromacs/bin/GMXRC
```


This will download and build first the prerequisite FFT library
followed by GROMACS. If you already have FFTW installed, you can remove that argument to "cmake". Overall, this build of GROMACS will be correct and reasonably fast on the machine upon which "cmake" ran.
On another machine, it may not run, or may not run fast. If you want to get the maximum value for your hardware with GROMACS, you will have to read further. Sadly, the interactions of hardware, libraries, and compilers are only going to continue to get more complex.

#### III. Examples
The example of CWBSol method is located in *‌cwbsol_example* folder. 
Running this example by command:
```
gmx_d cwbsol -solvent tip3p.gaff -solute ala2.solute -nr 30x30x30 -traj ala2.pdb -cmd closure=kh rism
```

## Deployment

CWBSol method is developed for Gromacs 2019 by [Prof. Xuhui Huang's group](http://compbio.ust.hk) 

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors
* **Prof. Xuhui Huang** - *Project leader* - [xuhuihuang](http://compbio.ust.hk/public_html/pmwiki-2.2.8/pmwiki.php?n=People.XuhuiHuang)
* **Dr. Siqin Cao** - *main developer* - [seechin](https://github.com/seechin)
* **Mr. Song Liu** - *developer* -[liusong299](https://github.com/liusong299/)

See also the list of [contributors](https://github.com/liusong299/gromacs-2019-CWBSol/graphs/contributors) who participated in this project.

## License

This project is licensed under the GPL License - see the [LICENSE](LICENSE) file for details

## Acknowledgments


## References:

1. Siqin Cao, Fu Kit Sheong, and Xuhui Huang, “Reference interaction site model with hydrophobicity induced density inhomogeneity: An analytical theory to compute solvation properties of large hydrophobic solutes in the mixture of polyatomic solvent molecules”, Journal of Chemical Physics 143, 054110 (2015)
2. Siqin Cao, Lizhe Zhu and Xuhui Huang, 3DRISM-HI-D2MSA: an improved analytic theory to compute solvent structure around hydrophobic solutes with proper treatment of solute–solvent electrostatic interactions, Mol. Phys. 116, 1003 (2017)
3. Siqin Cao, Kirill Konovalov, Ilona Christy Unarta and Xuhui Huang, Recent Developments in Integral Equation Theory for Solvation to Treat Density Inhomogeneity at Solute-Solvent Interface, Advanced Theory and Simulations (2019) DOI: 10.1002/adts.201900049