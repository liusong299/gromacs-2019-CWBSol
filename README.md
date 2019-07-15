# Gromacs CWBSol
The original release of CWBSol method with Gromacs 2019

## Authors: Cao Siqin (2015-2019)

## Getting Started
The 3D reference interaction site model (CWBSol) is a powerful tool to study the thermodynamic and structural properties of liquids. However, for hydrophobic solutes, the inhomogeneity of the solvent density around them poses a great challenge to the 3DRISM-HI theory. To address this issue, we have previously introduced the hydrophobic-induced density inhomogeneity theory (HI) for purely hydrophobic solutes. To further consider the complex hydrophobic solutes containing partial charges, here we propose the D2MSA closure to incorporate the short-range and long-range interactions with the D2 closure and the mean spherical approximation, respectively. We demonstrate that our new theory can compute the solvent distributions around real hydrophobic solutes in water and complex organic solvents that agree well with the explicit solvent molecular dynamics simulations.

### Prerequisites

#### I. Package requirement
* FFTW3 with double precision (required)
* GROMACS (optional, for reading XTC)
* pthread (optional, for pthread paralleling)
* LIBZ (optional, for output compression)

## References:

1. Siqin Cao, Fu Kit Sheong, and Xuhui Huang, “Reference interaction site model with hydrophobicity induced density inhomogeneity: An analytical theory to compute solvation properties of large hydrophobic solutes in the mixture of polyatomic solvent molecules”, Journal of Chemical Physics 143, 054110 (2015)
2. Siqin Cao, Lizhe Zhu and Xuhui Huang, 3DRISM-HI-D2MSA: an improved analytic theory to compute solvent structure around hydrophobic solutes with proper treatment of solute–solvent electrostatic interactions, Mol. Phys. 116, 1003 (2017)
3. Siqin Cao, Kirill Konovalov, Ilona Christy Unarta and Xuhui Huang, Recent Developments in Integral Equation Theory for Solvation to Treat Density Inhomogeneity at Solute-Solvent Interface, Advanced Theory and Simulations (2019) DOI: 10.1002/adts.201900049

# Project Title

One Paragraph of project description goes here

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
