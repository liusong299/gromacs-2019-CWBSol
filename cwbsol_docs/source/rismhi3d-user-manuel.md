<span>Version 0.243, (c) Siqin Cao</span>

Overview
========

The RISMHI3D is a software package of solvent structure calculation.
Based on the

is based on the integration equation theory (IET) <span>[cite]</span> of
liquid, the hydrophobicity induced density inhomogeneity (HI) theory
<span>[cite]</span> and

The RISMHI3D can perform the following calculations:

Theory
======

Installation guide
==================

A step-by-step guide of running RISMHI3D
========================================

Step 1. prepare the solvent setting file {#step-1.-prepare-the-solvent-setting-file .unnumbered}
----------------------------------------

To run the RISMHI3D to calculate the solvent distributions around given
solutes, three input files are required: the solvent setting file, the
solute forcefield file, and the solute conformations. Some solvent
setting file can be found in the software library; and if our library
doesn’t have your solvent, please generate a new solvent setting file
according to Sec. [sec:io.solvent].

Now suppose you have `tip3p.gaff` in `$IETLIB` folder.

Step 2. prepare the solute conformation file {#step-2.-prepare-the-solute-conformation-file .unnumbered}
--------------------------------------------

Then you need to find the file of solute conformation(s). This can be in
either PDB, GRO or XTC format (see Sec. [sec:io.traj] for details). The
conformation(s) should not contain any solvent of the kind that you want
to perform RISMHI3D for. The solute is not required to be electrically
neutralized.

Now suppose you have `rcsb.pdb` downloaded from RCSB, which may not be
compatible with RISMHI3D. There are two ways to change it to the
RISMHI3D compatible format. In the first way, you can do it youself by
keeping only the lines that begin with “ATOM”, and add one line to
define the box vector (suppose the box is
$100\times 110\times 120 {\rm\AA}^3$, and please ignore the rulers in
gray texts):

`|——–+———+———+———+———+———+———|`

`CRYST1  100.000  110.000  120.000  90.00  90.00  90.00 P 1           1`

`|——–+———+———+———+———+———+———|`

Step 3. prepare the solute forcefield file {#step-3.-prepare-the-solute-forcefield-file .unnumbered}
------------------------------------------

Then you need to prepare the solute forcefield file for this solute. The
solute forcefield file can be generated with the GROMACS top file (see
Sec. [sec:io.solute]), so before you can first generate the top file.
This can be done with either the GROMACS package or AMBER Tools and
ACPype. If you have GROMACS:

`gmx pdb2gmx -f rcsb.pdb -o solute.pdb -p solute.top`

or if you have both AMBER Tools and ACPype:

`antechamber -i solute.pdb -fi pdb -o solute.mol2 -fo mol2 -s 2 -c bcc`

`acpype.py -i solute.mol2 -c user`

`cp solute.acpype/solute_GMX.top solute.top`

`cp solute.acpype/solute_GMX.pdb solute.pdb`

Then, the solute forcefield file `soltue.ff` can be generated with the
following commands:

`GMXDATA=/opt/gromacs/share gmxtop2solute -top solute.top -o solute.ff`

or

`gmxtop2solute -top solute.top -ffpath /opt/gromacs/share/gromacs/top -o solute.ff`

(Please note that the atoms in the solute conformation file should be
exactly the same as in the solute forcefield file.)

Step 4. Run the 3DRIMSHI {#step-4.-run-the-3drimshi .unnumbered}
------------------------

A simple run of RISMHI3D is shown as below, the 3DRISM-KH calculation of
the TIP3P water around a solute molecule defined in output.ff and
output.pdb:

`rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -do-rism-kh -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60`

For advanced users: the above command is equivalent to the following
command:

`rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -cmd closure=kh rism report:energy display:rdf savee:cmd,guv -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60`

<span>(the green text are optional parameters and commands)</span>

![The screen output of running 3DRISMHI on alanine
dipeptide](rismhi3d-demo-screen.png "fig:") [fig:demo-screen]

The screen output can be seen in Fig. [fig:demo-screen]. At the same
time, a file `output.ts4s` containing the command, density distributions
and the HI density distributions is generated in the current folder.
This file can be decoded with `ts4sdump`:

`ts4sdump output.ts4s`

`ts4sdump -e 2 output.ts4s`

`ts4sdump -d 2 output.ts4s`

For hydrophobic solutes, the HI is suggested to be performed. For the
above example, the 3DRISM-HI-D2MSA calculation can be performed with the
following command:

`rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -do-rismhi-d2 -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60`

Again, for advanced users: this command is equivalent to the following
command:

`rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -cmd hshi closure=d2 rism report:energy display:rdf savee:cmd,guv -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60`

Input and output
================

RISMHI3D work with three input files and some options. The three input
files are:

-   Solvent settings (-p, -solvent, see Sec. [sec:io.solvent]): the
    information of te solvent

-   Solute forcefield (-s, -solute, see Sec. [sec:io.solute]):
    forcefield parameters of the solute

-   Solute conformation(s) (-f, -traj, see Sec. [sec:io.traj]): the
    atomic structure of solute molecules

-   Parameters, options and commands (see Sec. [sec:param]): commands to
    run, number of grids, number of threads, etc.

The output of RISMHI3D consists of screen outputs, TS4S outputs and RDF
outputs:

-   Screen outputs (see Sec. [sec:io.screen]): brief reports and debug
    information

-   TS4S files (see Sec. [sec:io.ts4s]): results in the grid space

-   RDF file (see Sec. [sec:io.rdf]): text file containing the
    user-defined RDF pairs

By default, only a few screen output will be generated. More details of
screen output, TS4S or RDF output will be shown by enabling options or
commands (see Sec. [sec:param])

Input: Solvent setting file {#sec:io.solvent}
---------------------------

The solvent setting file is specified with -p or -solvent. If the
current folder doesn’t have the solvent setting file, the program will
search for the path defined in \$IETLIB. Normally the filename should
not begin with “-”, and if you do have such a file, use “`–`” before the
filename (e.g. `-p – -solvent_settings.gaff`)

The solvent setting file contains several sections: [solvent], [atom],
[bond], and [gvv\_map] ([gvv\_map] is optional).

The <span>**[solvent]**</span> section mainly contains the following
parameters:

-   <span>**ff**</span>: the forcefield type specifier, can be:
    gaff/amber/opls. The mixing rule of VdW radiuses are arithmetic
    (-arith-sigma) in gaff/amber or geometric (-geo-sigma) in opls. The
    energy unit is kJ/mol (-Tdef 120.27) in gaff/opls or kcal/mol (-Tdef
    502.97) in amber.

-   <span>**rvdw, rcoul, rc**</span>: the cutoff distance (in nm) for
    the interactions. -rc will specify the cutoff for both -rcoul and
    -rvdw. The LJ interaction is calculated at a hard cutoff at -rvdw,
    while the Coulomb interaction can be calculated with either a PME
    swithed at -rcoul (if PME is enabled with -pme) or a hard cutoff at
    -rcoul (if PME is disabled by -nopme)

-   <span>**density**</span> and <span>**bulk-density**</span>: a list
    of the densities or bulk densities of all solvent components. The
    unit is nm$^{-3}$. The bulk densities are the pure liquid densities
    of each solvent component.

-   <span>**gvv**</span>: the solvent-solvent correlations.
    <span>**gvv**</span> is followed by the name of gvv file and the
    grid size, e.g. “-gvv 0.001 tip3p.gvv” or “-gvv tip3p.gvv 0.001”
    means that the gvv is defined in tip3p.gvv with the grid size of
    0.001 nm. The gvv file contains a number of columns, each column is
    the correlation function of certain pair of solvent sites. The
    mapping between columns and solvent site pairs are defined by
    [gvv\_map] section or in a default order (if [gvv\_map] section is
    missing). The real mapping of gvv can be seen when running the
    program with “-list” option.

-   <span>**zeta**</span>: the solvent-solvent zeta correlations, only
    used in HI calculations. The definition of zeta is similar to gvv,
    while the only difference is the unit of zeta is energy. zeta has no
    mapping between columns and molecule pairs, and should contain
    N$\times$N columns for N solvent molecules.

-   <span>**dielect**</span> or <span>**dipole**</span>: a list of
    dielectric constants (vacuum is 1) or dipole moments (unit is
    $\rm e\cdot nm$) of all molecules. The dielectric constants and
    dipole moments are optional and only used in certain algorithms
    (e.g. -Coulomb dielect, see Sec. [sec:param] for details)

The <span>**[atom]**</span> section defines all the atoms. The [atom]
section consists seven columns (more columns are ignored), separated by
spaces or tabs:

  --- ---------------
   1  atom name
   2  moleucle name
   3  index
   4  group index
   5  charge
   6  sigma (nm)
   7  epsilon
  --- ---------------

The group index is the site index, and the atoms with the same group (or
site) index can be treated as one site in RISM calculation. It is highly
suggested that completely equivalent atoms are grouped into one site,
which will greatly reduce the computational cost.

The <span>**[bond]**</span> section defines the bonds and pairs of each
pair of atoms. The <span>**[bond]**</span> section consists of three or
four columes, where the first two columns are the two atoms, and the
third column is the bond length or fixed pair distance between the pair
of atoms (unit: nm). The fourth column is the RMSD of the fluctuations
of bond lengths (or pair distances), which is optional and will be
treated as $0$ if the fourth column is missing. Warning: don’t define
one bond/pair twice.

Input: Solute forcefield {#sec:io.solute}
------------------------

The solute forcefield file is specified with -s or -solute. This file
contains only <span>**[solute]**</span> section(s). The
<span>**[solute]**</span> section contains six or eight (more column are
ignored) columns:

   column    6-col format     8-col format
  -------- ---------------- ----------------
     1        atom name        atom index
     2      moleucle name      atom name
     3        atom mass        residue id
     4      partial charge   moleucle name
     5        sigma (nm)       atom mass
     6         epsilon       partial charge
     7                         sigma (nm)
     8                          epsilon

The solute forcefield can be simply translated from the GROMACS top file
with the “gmxtop2solute” tool provided in the software package.

Input: Solute conformation(s) {#sec:io.traj}
-----------------------------

The solute conformation(s) are defined in the trajectory file, specified
with -f or -traj. The trajectory file can be a PDB, GRO or XTC file.

The PDB file: only lines begin with “ATOM” or “CRYST1” will be
processed. The number of atoms defined by the ATOM lines should be
consistent with the solute forcefield file. The box size should be
defined with a CRYST1 line. The PDB file can contain multiple frames,
separated by “ENDMDL”.

The GRO file can also contain multiple frames.

The XTC file can be processed only when the software is compiled with
\_GROMACS\_ options. It’s fine to turn off all \_GROMACS\_ flags when
compling the software, as long as you don’t have GROMACS or don’t want
to be bothered by this feature. The frames of XTC between the time (in
ps) defined in <span>**-b**</span> and <span>**-e**</span> are handled,
and <span>**-dt**</span> specifies the time interval between frames that
handled in RISMHI3D.

Output: screen report {#sec:io.screen}
---------------------

The screen output consists of the information of running and some brief
reports of calculation results.

The running information can be muted with <span>**-v 0**</span>. If you
want to monitor the status of running, please use <span>**-v 1**</span>,
<span>**-v 2**</span> or simply <span>**-v**</span>.

More running information is shown in debug mode. <span>**-debug
0**</span> will mute all debug message, and <span>**-debug 1**</span>
will allow to show some important messages (e.g. real location of input
files, allocated memory, thread, etc.). In <span>**-debug 2**</span>, a
detailed running process is printed on screen, including calling of
major functions in the source code. Further in <span>**-debug 3**</span>
or <span>**-debug-crc**</span>, the CRC check sum of important memory
bulks will be displayed on the screen. Please note that -debug 3 or
-debug-crc will perform the CRC calculation, which will require
additional computing time.

A detail report of time consumption will be displayed at the end with
<span>**-v 2**</span> or <span>**-debug 1/2/3**</span>.

The brief reports of calculation results will be displayed according to
the command defined in the command queue. Please see Sec. [sec:param]
for details of commands and the command queue. For example, if you want
a detailed report of energy or correlation functions, you need to add a
“-cmd report:energy” or “-cmd report:cuv” command; if you want to show
the rdf, you need to add a “-cmd display:rdf” command.

The screen output can be redirected to a log file specified with
<span>**-log**</span> (<span>**-log screen**</span>, <span>**-log
stdout**</span> and <span>**-log stderr**</span> will redirect the
screen output to stdout or stderr).

Output: TS4S file {#sec:io.ts4s}
-----------------

The TS4S file contains one or more frames of 4D tensors. Normally in
RISMHI3D, the 4D tensors can be LJ potentials, Coulomb potentials, total
correlation, direct correlation, total density profile and HI density.

The name of the output TS4S file is specified with
<span>**-o[v][0/1/2]**</span> or <span>**-a[0/1/2]**</span>, where “ov”
represepts overriding, and “a” represents appending. The TS4S data can
be uncompressed (<span>**-o[v]0**</span> or <span>**-a0**</span>) or
compressed with ZLIB (if you compiled with the \_ZLIB\_ option). Please
note that the TS4S file may be extrodinarily huge, and the
<span>**-significant-digits**</span> or <span>**-sd**</span> option may
be helpful for higher compressibility with lower accuracy of output data
(e.g. <span>**-sd 5**</span> will keep only five significant digits, and
<span>**-sd float**</span> will trim all output data to float).

The output of TS4S will be performed with the “save” command defined in
the command queue (see Sec. [sec:param] for details). The TS4S file will
be generated at the first saving command. The filename extension of TS4S
is always “.ts4s”. If the filename is not specified, a default filename
solute.solvent.YYMMDD\_HHMM.ts4s will be used.

Please use “ts4sdump” (provided in the software package) to check and
decode the TS4S file.

<span>*(For developer:) The file format of a TS4S file is defined in
compress.cpp of the source code. TS4S consists a number of data blocks,
each block corresponding to a 4D tensor. One 4D tensor block begins with
the IETSPageHeader structure, followed by comment text and tensor data.
The decoded 4D tensors are organized in the order of
“tensor[solvent\_site][z][y][x]”.*</span>

Output: RDF file {#sec:io.rdf}
----------------

The RDF can be calculated and displayed or written to a text file. The
RDF groups are defined in “<span>**-rdf-grps**</span> u1-v1,u1-v2,...”,
where u1, u2 ... stand for the indexes of solute atoms and v1, v2 ...
stand for the indexes of solvent sites. The number of RDF bins is
defined in “<span>**-rdf-bins**</span> 50”, where 50 is the default bin
number. The RDF will be performed up to the distance defined in
“<span>**-rc**</span>”. For an advanced user who wants to see the RDFs
of HI densities or direct correlations, the
“<span>**-rdf-content**</span>” option can be used to specify the RDF to
calculate: rdf (density profile, by default), h (direct correlation), dd
(HI density), c (direct correlation), ch (ch=$c*h$), lj (LJ potential)
or coul (Coulomb potential).

The displaying or saving of RDF is performed in “-cmd display:rdf” or
“-cmd save:rdf”. The RDF is calculated when necessary, so no RDF
calculation will be performed if the RDF is not displayed or saved.
“-cmd display:rdf” will display the RDF on screen (or -log file), while
“-cmd save:rdf” will save the RDF to a text file. The name of the RDF
file begins with the name of TS4S file, and ended with a “.rdf”
extension.

Parameters, options and commands {#sec:param}
================================

Parameters for running {#sec:param.param}
----------------------

In the command line parameters, you can specify the running parameters,
tell the software what to compute and display, and override some solvent
settings. (Please note that all the settings in the
<span>**[solvent]**</span> section can be overridden by the command line
parameters.) Here is a list of the major parameters for running:

-   <span>**-nt**</span> and <span>**-np**</span>: number of parallel
    runnings. -np 1 or -nt 1 will disable the paralleling.
    <span>**-nt**</span> will use multithread while <span>**-np**</span>
    will use fork. Please note that paralleling is still possible (with
    -np) even if you hate or don’t have pthread. <span>**-np**</span> is
    enabled with \_LOCALPARALLEL\_ option in compiling, while
    <span>**-nt**</span> requires both \_LOCALPARALLEL\_ and
    \_LOCALPARALLEL\_PTHREAD\_ in compiling. Note: the RISMHI3D will
    automatically set -nt to the maximum number of your CPU cores. So
    don’t forget to change -nt or -np if you are running RISMHI3D on a
    cluster.

-   <span>**-nice**</span>: the nice level of running this software.

-   <span>**-nr **</span>: the grid number. It can be one number for
    both three dimensions; three numbers for X, Y and Z respectively; or
    NXxNYxNZ. e.g. “-nr 50” is equivalent to “-nr 50 50 50” and “-nr
    50x50x50”. The grids of X, Y and Z can be different, e.g. “-nr
    50x60x70”. The grid numbers are suggested to be even numbers.

-   <span>**-step**</span> 100: set the maximum steps for HI and RISM.

-   <span>**-interact**</span>: enable the interactive mode. Pressing
    enter at running will halt the calculation, and you can print the
    report, change parameters, continue running, end current RISM/HI
    calculation or quite the program. Warning: don’t enable the
    interactive mode when you are running RISMHI3D in background (e.g.
    nohup or on a cluster).

-   <span>**-do-rism-kh**</span>, <span>**-do-rismhi-kh**</span> and
    <span>**-do-rismhi-d2**</span>: perform 3DRISM or 3DRISMHI with the
    KH or D2MSA closures. The RDFs will be displayed on screen if
    <span>**-rdf-grps**</span> are defined, while a TS4S file containing
    the running command and the density profile will be generated if the
    output file is specified.

-   <span>**-cmd**</span> or <span>**-do**</span>: the command(s) you
    want to run. The software will do nothing if the command queue is
    empty. See Sec. [sec:param.command] for details.

And here is a list of some important <span>**[solvent]**</span> settings
that can be redefined:

-   <span>**-temperature**</span> 298: specify the temperature in Kelvin

-   <span>**-Tdef**</span> 120.27: specify the temperature of defining
    the energy unit. 1 kJ/mol = 120.27 K, 1 kcal/mol = 502.97 K. Can be
    specified with <span>**-ff**</span>.

-   <span>**-arith-sigma**</span> and <span>**-geo-sigma**</span>: set
    the combing rule of VdW raiuds to arithmetic or geometric averaging.

-   <span>**-density**</span> and <span>**-bulk-density**</span> can be
    redefined in the command line parameters.

-   <span>**-ndiis**</span>, <span>**-ndiisrism**</span> and
    <span>**-ndiishi**</span>: the maximum steps of DIIS. More steps of
    DIIS will have better convergence of self-consistent-field
    iterations while require more memory.

-   <span>**-delvv**</span>, <span>**-delrism**</span> and
    <span>**-delhi**</span>: the step in factor of the
    self-consistent-field iterations. <span>**-delvv**</span> for both
    RISM and HI. Both are 1 by default, and 0.7 is recommended by RISM
    <span>[cite]</span>.

-   <span>**-errtol**</span>, <span>**-errtolrism**</span> and
    <span>**-errtolhi**</span>: the error tolerance of convergence. Both
    are $10^{-12}$ by default. Although this is fine for HI, the RISM
    iterations have little chances to reach the error of $10^{-12}$.
    This number is suggested to be $10^{-7}$ in AMBER RISM.

-   <span>**-bound-to-ram**</span> or
    <span>**-ignore-memory-capacity**</span>: by default, the software
    will detect the capacity of the physical memory, and will terminate
    when the memory is exceeding the physical memory. This feature can
    be disabled by <span>**-ignore-memory-capacity**</span>, which will
    cause extremely low computational efficiency as well as high risk to
    damage your hard disk. Ignore the memory capacity check only when
    you know what you are doing.

Other advanced options can be seen in Sec. [sec:param.options]

Commands and command queue {#sec:param.command}
--------------------------

The command forms a command queue, and will be performed one by one
after the frames are read from the trajectory file. The command(s)
specified in <span>**-cmd**</span> will be added to the command queue.
You can use to relocate the current command instead of inserting it to
the end of the queue, e.g. “-cmd report@5:energy” will insert “-cmd
report:energy” to the 5th command (the command queue begins with 1). b
and e can be used to force the command to runs before handling of any
frame or after handling all frames, e.g. “-cmd report:rdf” will report
the RDF at each frame, and “-cmd report@end:rdf” will report the overall
RDF after all frames have been processed.

Followings are some basic commands for setting closures, running HI/RISM
calculations, and generating reports or output files of the output
results:

-   <span>**build-ff**</span>: <span>**-cmd build-ff**</span>: force to
    rebuild the forcefield. The forcefield is built automatically with
    <span>**rism**</span> or <span>hshi</span>, and this command is used
    only when you need to do something without perfoming RISM or HI.

-   <span>**rism**</span>: <span>**-cmd rism**</span>,step=100: perform
    3DRISM with specified number of maximum steps and the closure(s)
    defined in the <span>**closure**</span> command before.

-   <span>**ssoz**</span>: <span>**-cmd ssoz**</span>,step=100: perform
    unrenormalized 3DRISM with specified number of maximum steps and the
    closure(s) defined in the <span>**closure**</span> command before.
    The unrenormalized 3DRISM is the very original version of 3DRISM,
    which has big issues in electrostatic interactions.

-   <span>**hshi**</span>: <span>**-cmd hshi**</span>,step=100: perform
    HI with specified number of maximum steps.

-   <span>**closure**</span>: <span>**-cmd
    closure=**</span>closure\_A[,closure\_B,...]: set the closure(s) for
    each molecule. Different sites can be calculated with different
    closures, and at most 20 closures are allowed to be specified here.
    Particularly, all the sovlent sites will use the same specified
    closure if only one closure is given here.

-   <span>**closure-a**</span>: set the closure(s) for each site instead
    of molecule.

-   <span>**closure-factor**</span>: <span>**-cmd
    cf=**</span>cf\_A[,cf\_B,...]: set the extra parameters that used in
    closures.

-   <span>**report**</span>: <span>**-cmd
    report:**</span>Euv/energy/cuv/rdf: generate a report on screen (or
    to -log file). <span>**-cmd report:Euv**</span> will display a brief
    report of total energies, while <span>**-cmd report:energy**</span>
    will display a detailed report of total energies. <span>**-cmd
    report:cuv**</span> will display the total direct correlations, and
    <span>**-cmd report:rdf**</span> will display the RDFs if you have
    defined the RDF groups. In addition, <span>**-cmd
    report:energy,cuv**</span> is equivalent to <span>**-cmd
    report:all**</span>.

-   <span>**display**</span>: <span>**-cmd
    display:**</span>lj/coul/Euv/energy/cuv/dN/dN0/TS/GGF/rdf: display
    the values of the chosen variables

-   <span>**save**</span>: <span>**-cmd
    save:**</span>cmd/lj/coul/cuv/huv/hlr/dd/ddp/nphi/guv/rmin/rdf: save
    the specified quantity. <span>**-cmd save:cmd**</span> will save the
    command line arguments to the TS4S file, and <span>**-cmd
    save:guv**</span> will save the density distributions of each
    solvent site at each spacial grid to the TS4S file. <span>**-cmd
    save:lj,coul,cuv,huv,hlr,dd**</span> will save LJ potentials,
    Coulomb potentials, direct correlations, total correlations, long
    range total correltations, HI density to the TS4S file, and
    <span>**-cmd save:rmin**</span> will save the minimal-to-solute
    distances to the TS4S file. Additionally, <span>**-cmd
    save:rdf**</span> will save the RDF to the RDF file (which is a text
    file).

-   <span>**savee**</span>: mostly the save as <span>**save**</span>.
    The only difference is that <span>**savee**</span> will perform
    saving only when the output TS4S file is explicitly specified, while
    <span>**save**</span> will use an default filename if the output
    TS4S file is not explicitly specified.

Addtionally, some shortcuts can be used to add a bundle of commands.
Don’t use two or more shortcuts, otherwise both the two sets of commands
will be performed.

-   <span>**-do-rism-kh**</span> = -cmd closure=kh rism report:energy
    display:rdf savee:cmd,guv

-   <span>**-do-rismhi-d2**</span> = -cmd hshi closure=d2 rism
    report:energy display:rdf savee:cmd,guv

-   <span>**-do-rismhi-kh**</span> = -cmd hshi closure=kh rism
    report:energy display:rdf savee:cmd,guv

Other options {#sec:param.options}
-------------

Below is a list of major advanced options that can be defined in both
the <span>**[solvent]**</span> section of the solvent setting file or
the command line parameters:

-   <span>**-lsa**</span> and <span>**-lsb**</span> for HI: the two
    parameters, $A$ and $B$ of the liquid equation of state
    <span>[cite]</span>. $B$ is automatically computed, and $A$ can be
    defined with <span>-lse-a</span> or <span>**-lsa**</span>. The
    recommended values of $A$ can be seen from previous experimental
    measurements <span>[cite]</span>.

-   <span>**-theta**</span> for HI: define the energy cutoff (in kT),
    above which the regions are treated as hard cores or no solvent
    regions. The default cutoff is 5 kT.

-   <span>**-Coulomb**</span>: the preprocessing algorithm of the
    Coulomb interactions in HI and RISM. Can be: none (=-Coulomb),
    dielect, or YukawaFFT.

-   <span>**-Yukawa**</span> 0.5: the same as <span>**-Coulomb
    YukawaFFT**</span> 0.5. 0.5 here is the characteristic length (unit:
    nm) of the exponential function of the Yukawa potential. The Debye
    length of the homogeneous liquid will be used if the characteristic
    length is not specified. The dielectric constant for the Yukawa
    potential is defined in <span>**-dielect-y**</span>.

-   <span>**-ccutoff**</span>: the cutoff of closures, which is used in
    the closures like PLHNC.

-   <span>**-rb**</span> or <span>**-Bohr-radius**</span>: the minimal
    raiuds of an atom. By default it is 0.052911 (nm).

-   <span>**-sd**</span> or <span>**-significant-digits**</span>: the
    significant digits in the TS4S files, can be “float”, “double” or a
    number (significant digits in decimal)

-   <span>**-closure-enhance-level**</span> 1 or
    <span>**-enhance-closure**</span> 1: scale down the changes of SCF
    iterations with $(1+h^2)^{\alpha/2}$, where $\alpha$ is the closure
    enhancement level. This option greatly helps the convergence of the
    self-consistent-field iterations and is turned on by default. Use
    <span>**-no-closure-enhance**</span> to turn this feature off.

-   <span>**-bounded-to-ram**</span>: don’t exceed the physical memory
    capacity when using memory. This feature is on by default, and use
    <span>**-ignore-memory-capacity**</span> to turn this feature off.

-   <span>**-xvv-extend**</span> 0: extend the <span>**gvv**</span> of
    solvent. If the input gvv contains the RDFs of $2$ nm, then
    -xvv-extend 5 will extend it to $10$ nm by filling the extended
    regions with $1$. This option is helpful when you
    <span>**gvv**</span> is poor, but don’t expect too much.


