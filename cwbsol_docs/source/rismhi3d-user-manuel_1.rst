Version 0.243, (c) Siqin Cao

Overview
========

The RISMHI3D is a software package of solvent structure calculation.
Based on the

is based on the integration equation theory (IET) [cite] of liquid, the
hydrophobicity induced density inhomogeneity (HI) theory [cite] and

The RISMHI3D can perform the following calculations:

Theory
======

Installation guide
==================

A step-by-step guide of running RISMHI3D
========================================

Step 1. prepare the solvent setting file
----------------------------------------

To run the RISMHI3D to calculate the solvent distributions around given
solutes, three input files are required: the solvent setting file, the
solute forcefield file, and the solute conformations. Some solvent
setting file can be found in the software library; and if our library
doesn’t have your solvent, please generate a new solvent setting file
according to Sec. [sec:io.solvent].

Now suppose you have ``tip3p.gaff`` in ``$IETLIB`` folder.

Step 2. prepare the solute conformation file
--------------------------------------------

Then you need to find the file of solute conformation(s). This can be in
either PDB, GRO or XTC format (see Sec. [sec:io.traj] for details). The
conformation(s) should not contain any solvent of the kind that you want
to perform RISMHI3D for. The solute is not required to be electrically
neutralized.

Now suppose you have ``rcsb.pdb`` downloaded from RCSB, which may not be
compatible with RISMHI3D. There are two ways to change it to the
RISMHI3D compatible format. In the first way, you can do it youself by
keeping only the lines that begin with “ATOM”, and add one line to
define the box vector (suppose the box is
:math:`100\times 110\times 120 {\rm\AA}^3`, and please ignore the rulers
in gray texts):

``|——–+———+———+———+———+———+———|``

``CRYST1  100.000  110.000  120.000  90.00  90.00  90.00 P 1           1``

``|——–+———+———+———+———+———+———|``

Step 3. prepare the solute forcefield file
------------------------------------------

Then you need to prepare the solute forcefield file for this solute. The
solute forcefield file can be generated with the GROMACS top file (see
Sec. [sec:io.solute]), so before you can first generate the top file.
This can be done with either the GROMACS package or AMBER Tools and
ACPype. If you have GROMACS:

``gmx pdb2gmx -f rcsb.pdb -o solute.pdb -p solute.top``

or if you have both AMBER Tools and ACPype:

``antechamber -i solute.pdb -fi pdb -o solute.mol2 -fo mol2 -s 2 -c bcc``

``acpype.py -i solute.mol2 -c user``

``cp solute.acpype/solute_GMX.top solute.top``

``cp solute.acpype/solute_GMX.pdb solute.pdb``

Then, the solute forcefield file ``soltue.ff`` can be generated with the
following commands:

``GMXDATA=/opt/gromacs/share gmxtop2solute -top solute.top -o solute.ff``

or

``gmxtop2solute -top solute.top -ffpath /opt/gromacs/share/gromacs/top -o solute.ff``

(Please note that the atoms in the solute conformation file should be
exactly the same as in the solute forcefield file.)

Step 4. Run the 3DRIMSHI
------------------------

A simple run of RISMHI3D is shown as below, the 3DRISM-KH calculation of
the TIP3P water around a solute molecule defined in output.ff and
output.pdb:

``rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -do-rism-kh -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60``

For advanced users: the above command is equivalent to the following
command:

``rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -cmd closure=kh rism report:energy display:rdf savee:cmd,guv -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60``

(the green text are optional parameters and commands)

|The screen output of running 3DRISMHI on alanine dipeptide|
[fig:demo-screen]

The screen output can be seen in Fig. [fig:demo-screen]. At the same
time, a file ``output.ts4s`` containing the command, density
distributions and the HI density distributions is generated in the
current folder. This file can be decoded with ``ts4sdump``:

``ts4sdump output.ts4s``

``ts4sdump -e 2 output.ts4s``

``ts4sdump -d 2 output.ts4s``

For hydrophobic solutes, the HI is suggested to be performed. For the
above example, the 3DRISM-HI-D2MSA calculation can be performed with the
following command:

``rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -do-rismhi-d2 -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60``

Again, for advanced users: this command is equivalent to the following
command:

``rismhi3d -p $IETLIB/tip3p.gaff -s solute.ff -f solute.pdb -nr 100x110x120 -rc 1.2 -cmd hshi closure=d2 rism report:energy display:rdf savee:cmd,guv -step 200 -o output -rdf-grps 2-1,2-2,9-1,9-2,19-1,19-2 -rdf-bins 60``

Input and output
================

RISMHI3D work with three input files and some options. The three input
files are:

-  Solvent settings (-p, -solvent, see Sec. [sec:io.solvent]): the
   information of te solvent

-  Solute forcefield (-s, -solute, see Sec. [sec:io.solute]): forcefield
   parameters of the solute

-  Solute conformation(s) (-f, -traj, see Sec. [sec:io.traj]): the
   atomic structure of solute molecules

-  Parameters, options and commands (see Sec. [sec:param]): commands to
   run, number of grids, number of threads, etc.

The output of RISMHI3D consists of screen outputs, TS4S outputs and RDF
outputs:

-  Screen outputs (see Sec. [sec:io.screen]): brief reports and debug
   information

-  TS4S files (see Sec. [sec:io.ts4s]): results in the grid space

-  RDF file (see Sec. [sec:io.rdf]): text file containing the
   user-defined RDF pairs

By default, only a few screen output will be generated. More details of
screen output, TS4S or RDF output will be shown by enabling options or
commands (see Sec. [sec:param])

Input: Solvent setting file
---------------------------

The solvent setting file is specified with -p or -solvent. If the
current folder doesn’t have the solvent setting file, the program will
search for the path defined in $IETLIB. Normally the filename should not
begin with “-”, and if you do have such a file, use “\ ``–``\ ” before
the filename (e.g. ``-p – -solvent_settings.gaff``)

The solvent setting file contains several sections: [solvent], [atom],
[bond], and [gvv\_map] ([gvv\_map] is optional).

The **[solvent]** section mainly contains the following parameters:

-  **ff**: the forcefield type specifier, can be: gaff/amber/opls. The
   mixing rule of VdW radiuses are arithmetic (-arith-sigma) in
   gaff/amber or geometric (-geo-sigma) in opls. The energy unit is
   kJ/mol (-Tdef 120.27) in gaff/opls or kcal/mol (-Tdef 502.97) in
   amber.

-  **rvdw, rcoul, rc**: the cutoff distance (in nm) for the
   interactions. -rc will specify the cutoff for both -rcoul and -rvdw.
   The LJ interaction is calculated at a hard cutoff at -rvdw, while the
   Coulomb interaction can be calculated with either a PME swithed at
   -rcoul (if PME is enabled with -pme) or a hard cutoff at -rcoul (if
   PME is disabled by -nopme)

-  **density** and **bulk-density**: a list of the densities or bulk
   densities of all solvent components. The unit is nm\ :math:`^{-3}`.
   The bulk densities are the pure liquid densities of each solvent
   component.

-  **gvv**: the solvent-solvent correlations. **gvv** is followed by the
   name of gvv file and the grid size, e.g. “-gvv 0.001 tip3p.gvv” or
   “-gvv tip3p.gvv 0.001” means that the gvv is defined in tip3p.gvv
   with the grid size of 0.001 nm. The gvv file contains a number of
   columns, each column is the correlation function of certain pair of
   solvent sites. The mapping between columns and solvent site pairs are
   defined by [gvv\_map] section or in a default order (if [gvv\_map]
   section is missing). The real mapping of gvv can be seen when running
   the program with “-list” option.

-  **zeta**: the solvent-solvent zeta correlations, only used in HI
   calculations. The definition of zeta is similar to gvv, while the
   only difference is the unit of zeta is energy. zeta has no mapping
   between columns and molecule pairs, and should contain
   N\ :math:`\times`\ N columns for N solvent molecules.

-  **dielect** or **dipole**: a list of dielectric constants (vacuum is
   1) or dipole moments (unit is :math:`\rm e\cdot nm`) of all
   molecules. The dielectric constants and dipole moments are optional
   and only used in certain algorithms (e.g. -Coulomb dielect, see Sec.
   [sec:param] for details)

The **[atom]** section defines all the atoms. The [atom] section
consists seven columns (more columns are ignored), separated by spaces
or tabs:

+-----+-----------------+
| 1   | atom name       |
+-----+-----------------+
| 2   | moleucle name   |
+-----+-----------------+
| 3   | index           |
+-----+-----------------+
| 4   | group index     |
+-----+-----------------+
| 5   | charge          |
+-----+-----------------+
| 6   | sigma (nm)      |
+-----+-----------------+
| 7   | epsilon         |
+-----+-----------------+

The group index is the site index, and the atoms with the same group (or
site) index can be treated as one site in RISM calculation. It is highly
suggested that completely equivalent atoms are grouped into one site,
which will greatly reduce the computational cost.

The **[bond]** section defines the bonds and pairs of each pair of
atoms. The **[bond]** section consists of three or four columes, where
the first two columns are the two atoms, and the third column is the
bond length or fixed pair distance between the pair of atoms (unit: nm).
The fourth column is the RMSD of the fluctuations of bond lengths (or
pair distances), which is optional and will be treated as :math:`0` if
the fourth column is missing. Warning: don’t define one bond/pair twice.

Input: Solute forcefield
------------------------

The solute forcefield file is specified with -s or -solute. This file
contains only **[solute]** section(s). The **[solute]** section contains
six or eight (more column are ignored) columns:

+----------+------------------+------------------+
| column   | 6-col format     | 8-col format     |
+==========+==================+==================+
| 1        | atom name        | atom index       |
+----------+------------------+------------------+
| 2        | moleucle name    | atom name        |
+----------+------------------+------------------+
| 3        | atom mass        | residue id       |
+----------+------------------+------------------+
| 4        | partial charge   | moleucle name    |
+----------+------------------+------------------+
| 5        | sigma (nm)       | atom mass        |
+----------+------------------+------------------+
| 6        | epsilon          | partial charge   |
+----------+------------------+------------------+
| 7        |                  | sigma (nm)       |
+----------+------------------+------------------+
| 8        |                  | epsilon          |
+----------+------------------+------------------+

The solute forcefield can be simply translated from the GROMACS top file
with the “gmxtop2solute” tool provided in the software package.

Input: Solute conformation(s)
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
ps) defined in **-b** and **-e** are handled, and **-dt** specifies the
time interval between frames that handled in RISMHI3D.

Output: screen report
---------------------

The screen output consists of the information of running and some brief
reports of calculation results.

The running information can be muted with **-v 0**. If you want to
monitor the status of running, please use **-v 1**, **-v 2** or simply
**-v**.

More running information is shown in debug mode. **-debug 0** will mute
all debug message, and **-debug 1** will allow to show some important
messages (e.g. real location of input files, allocated memory, thread,
etc.). In **-debug 2**, a detailed running process is printed on screen,
including calling of major functions in the source code. Further in
**-debug 3** or **-debug-crc**, the CRC check sum of important memory
bulks will be displayed on the screen. Please note that -debug 3 or
-debug-crc will perform the CRC calculation, which will require
additional computing time.

A detail report of time consumption will be displayed at the end with
**-v 2** or **-debug 1/2/3**.

The brief reports of calculation results will be displayed according to
the command defined in the command queue. Please see Sec. [sec:param]
for details of commands and the command queue. For example, if you want
a detailed report of energy or correlation functions, you need to add a
“-cmd report:energy” or “-cmd report:cuv” command; if you want to show
the rdf, you need to add a “-cmd display:rdf” command.

The screen output can be redirected to a log file specified with
**-log** (**-log screen**, **-log stdout** and **-log stderr** will
redirect the screen output to stdout or stderr).

Output: TS4S file
-----------------

The TS4S file contains one or more frames of 4D tensors. Normally in
RISMHI3D, the 4D tensors can be LJ potentials, Coulomb potentials, total
correlation, direct correlation, total density profile and HI density.

The name of the output TS4S file is specified with **-o[v][0/1/2]** or
**-a[0/1/2]**, where “ov” represepts overriding, and “a” represents
appending. The TS4S data can be uncompressed (**-o[v]0** or **-a0**) or
compressed with ZLIB (if you compiled with the \_ZLIB\_ option). Please
note that the TS4S file may be extrodinarily huge, and the
**-significant-digits** or **-sd** option may be helpful for higher
compressibility with lower accuracy of output data (e.g. **-sd 5** will
keep only five significant digits, and **-sd float** will trim all
output data to float).

The output of TS4S will be performed with the “save” command defined in
the command queue (see Sec. [sec:param] for details). The TS4S file will
be generated at the first saving command. The filename extension of TS4S
is always “.ts4s”. If the filename is not specified, a default filename
solute.solvent.YYMMDD\_HHMM.ts4s will be used.

Please use “ts4sdump” (provided in the software package) to check and
decode the TS4S file.

*(For developer:) The file format of a TS4S file is defined in
compress.cpp of the source code. TS4S consists a number of data blocks,
each block corresponding to a 4D tensor. One 4D tensor block begins with
the IETSPageHeader structure, followed by comment text and tensor data.
The decoded 4D tensors are organized in the order of
“tensor[solvent\_site][z][y][x]”.*

Output: RDF file
----------------

The RDF can be calculated and displayed or written to a text file. The
RDF groups are defined in “**-rdf-grps** u1-v1,u1-v2,...”, where u1, u2
... stand for the indexes of solute atoms and v1, v2 ... stand for the
indexes of solvent sites. The number of RDF bins is defined in
“**-rdf-bins** 50”, where 50 is the default bin number. The RDF will be
performed up to the distance defined in “**-rc**”. For an advanced user
who wants to see the RDFs of HI densities or direct correlations, the
“**-rdf-content**” option can be used to specify the RDF to calculate:
rdf (density profile, by default), h (direct correlation), dd (HI
density), c (direct correlation), ch (ch=:math:`c*h`), lj (LJ potential)
or coul (Coulomb potential).

The displaying or saving of RDF is performed in “-cmd display:rdf” or
“-cmd save:rdf”. The RDF is calculated when necessary, so no RDF
calculation will be performed if the RDF is not displayed or saved.
“-cmd display:rdf” will display the RDF on screen (or -log file), while
“-cmd save:rdf” will save the RDF to a text file. The name of the RDF
file begins with the name of TS4S file, and ended with a “.rdf”
extension.

Parameters, options and commands
================================

Parameters for running
----------------------

In the command line parameters, you can specify the running parameters,
tell the software what to compute and display, and override some solvent
settings. (Please note that all the settings in the **[solvent]**
section can be overridden by the command line parameters.) Here is a
list of the major parameters for running:

-  **-nt** and **-np**: number of parallel runnings. -np 1 or -nt 1 will
   disable the paralleling. **-nt** will use multithread while **-np**
   will use fork. Please note that paralleling is still possible (with
   -np) even if you hate or don’t have pthread. **-np** is enabled with
   \_LOCALPARALLEL\_ option in compiling, while **-nt** requires both
   \_LOCALPARALLEL\_ and \_LOCALPARALLEL\_PTHREAD\_ in compiling. Note:
   the RISMHI3D will automatically set -nt to the maximum number of your
   CPU cores. So don’t forget to change -nt or -np if you are running
   RISMHI3D on a cluster.

-  **-nice**: the nice level of running this software.

-  **-nr **: the grid number. It can be one number for both three
   dimensions; three numbers for X, Y and Z respectively; or NXxNYxNZ.
   e.g. “-nr 50” is equivalent to “-nr 50 50 50” and “-nr 50x50x50”. The
   grids of X, Y and Z can be different, e.g. “-nr 50x60x70”. The grid
   numbers are suggested to be even numbers.

-  **-step** 100: set the maximum steps for HI and RISM.

-  **-interact**: enable the interactive mode. Pressing enter at running
   will halt the calculation, and you can print the report, change
   parameters, continue running, end current RISM/HI calculation or
   quite the program. Warning: don’t enable the interactive mode when
   you are running RISMHI3D in background (e.g. nohup or on a cluster).

-  **-do-rism-kh**, **-do-rismhi-kh** and **-do-rismhi-d2**: perform
   3DRISM or 3DRISMHI with the KH or D2MSA closures. The RDFs will be
   displayed on screen if **-rdf-grps** are defined, while a TS4S file
   containing the running command and the density profile will be
   generated if the output file is specified.

-  **-cmd** or **-do**: the command(s) you want to run. The software
   will do nothing if the command queue is empty. See Sec.
   [sec:param.command] for details.

And here is a list of some important **[solvent]** settings that can be
redefined:

-  **-temperature** 298: specify the temperature in Kelvin

-  **-Tdef** 120.27: specify the temperature of defining the energy
   unit. 1 kJ/mol = 120.27 K, 1 kcal/mol = 502.97 K. Can be specified
   with **-ff**.

-  **-arith-sigma** and **-geo-sigma**: set the combing rule of VdW
   raiuds to arithmetic or geometric averaging.

-  **-density** and **-bulk-density** can be redefined in the command
   line parameters.

-  **-ndiis**, **-ndiisrism** and **-ndiishi**: the maximum steps of
   DIIS. More steps of DIIS will have better convergence of
   self-consistent-field iterations while require more memory.

-  **-delvv**, **-delrism** and **-delhi**: the step in factor of the
   self-consistent-field iterations. **-delvv** for both RISM and HI.
   Both are 1 by default, and 0.7 is recommended by RISM [cite].

-  **-errtol**, **-errtolrism** and **-errtolhi**: the error tolerance
   of convergence. Both are :math:`10^{-12}` by default. Although this
   is fine for HI, the RISM iterations have little chances to reach the
   error of :math:`10^{-12}`. This number is suggested to be
   :math:`10^{-7}` in AMBER RISM.

-  **-bound-to-ram** or **-ignore-memory-capacity**: by default, the
   software will detect the capacity of the physical memory, and will
   terminate when the memory is exceeding the physical memory. This
   feature can be disabled by **-ignore-memory-capacity**, which will
   cause extremely low computational efficiency as well as high risk to
   damage your hard disk. Ignore the memory capacity check only when you
   know what you are doing.

Other advanced options can be seen in Sec. [sec:param.options]

Commands and command queue
--------------------------

The command forms a command queue, and will be performed one by one
after the frames are read from the trajectory file. The command(s)
specified in **-cmd** will be added to the command queue. You can use to
relocate the current command instead of inserting it to the end of the
queue, e.g. “-cmd report@5:energy” will insert “-cmd report:energy” to
the 5th command (the command queue begins with 1). b and e can be used
to force the command to runs before handling of any frame or after
handling all frames, e.g. “-cmd report:rdf” will report the RDF at each
frame, and “-cmd report@end:rdf” will report the overall RDF after all
frames have been processed.

Followings are some basic commands for setting closures, running HI/RISM
calculations, and generating reports or output files of the output
results:

-  **build-ff**: **-cmd build-ff**: force to rebuild the forcefield. The
   forcefield is built automatically with **rism** or hshi, and this
   command is used only when you need to do something without perfoming
   RISM or HI.

-  **rism**: **-cmd rism**,step=100: perform 3DRISM with specified
   number of maximum steps and the closure(s) defined in the **closure**
   command before.

-  **ssoz**: **-cmd ssoz**,step=100: perform unrenormalized 3DRISM with
   specified number of maximum steps and the closure(s) defined in the
   **closure** command before. The unrenormalized 3DRISM is the very
   original version of 3DRISM, which has big issues in electrostatic
   interactions.

-  **hshi**: **-cmd hshi**,step=100: perform HI with specified number of
   maximum steps.

-  **closure**: **-cmd closure=**closure\_A[,closure\_B,...]: set the
   closure(s) for each molecule. Different sites can be calculated with
   different closures, and at most 20 closures are allowed to be
   specified here. Particularly, all the sovlent sites will use the same
   specified closure if only one closure is given here.

-  **closure-a**: set the closure(s) for each site instead of molecule.

-  **closure-factor**: **-cmd cf=**cf\_A[,cf\_B,...]: set the extra
   parameters that used in closures.

-  **report**: **-cmd report:**Euv/energy/cuv/rdf: generate a report on
   screen (or to -log file). **-cmd report:Euv** will display a brief
   report of total energies, while **-cmd report:energy** will display a
   detailed report of total energies. **-cmd report:cuv** will display
   the total direct correlations, and **-cmd report:rdf** will display
   the RDFs if you have defined the RDF groups. In addition, **-cmd
   report:energy,cuv** is equivalent to **-cmd report:all**.

-  **display**: **-cmd
   display:**lj/coul/Euv/energy/cuv/dN/dN0/TS/GGF/rdf: display the
   values of the chosen variables

-  **save**: **-cmd
   save:**cmd/lj/coul/cuv/huv/hlr/dd/ddp/nphi/guv/rmin/rdf: save the
   specified quantity. **-cmd save:cmd** will save the command line
   arguments to the TS4S file, and **-cmd save:guv** will save the
   density distributions of each solvent site at each spacial grid to
   the TS4S file. **-cmd save:lj,coul,cuv,huv,hlr,dd** will save LJ
   potentials, Coulomb potentials, direct correlations, total
   correlations, long range total correltations, HI density to the TS4S
   file, and **-cmd save:rmin** will save the minimal-to-solute
   distances to the TS4S file. Additionally, **-cmd save:rdf** will save
   the RDF to the RDF file (which is a text file).

-  **savee**: mostly the save as **save**. The only difference is that
   **savee** will perform saving only when the output TS4S file is
   explicitly specified, while **save** will use an default filename if
   the output TS4S file is not explicitly specified.

Addtionally, some shortcuts can be used to add a bundle of commands.
Don’t use two or more shortcuts, otherwise both the two sets of commands
will be performed.

-  **-do-rism-kh** = -cmd closure=kh rism report:energy display:rdf
   savee:cmd,guv

-  **-do-rismhi-d2** = -cmd hshi closure=d2 rism report:energy
   display:rdf savee:cmd,guv

-  **-do-rismhi-kh** = -cmd hshi closure=kh rism report:energy
   display:rdf savee:cmd,guv

Other options
-------------

Below is a list of major advanced options that can be defined in both
the **[solvent]** section of the solvent setting file or the command
line parameters:

-  **-lsa** and **-lsb** for HI: the two parameters, :math:`A` and
   :math:`B` of the liquid equation of state [cite]. :math:`B` is
   automatically computed, and :math:`A` can be defined with -lse-a or
   **-lsa**. The recommended values of :math:`A` can be seen from
   previous experimental measurements [cite].

-  **-theta** for HI: define the energy cutoff (in kT), above which the
   regions are treated as hard cores or no solvent regions. The default
   cutoff is 5 kT.

-  **-Coulomb**: the preprocessing algorithm of the Coulomb interactions
   in HI and RISM. Can be: none (=-Coulomb), dielect, or YukawaFFT.

-  **-Yukawa** 0.5: the same as **-Coulomb YukawaFFT** 0.5. 0.5 here is
   the characteristic length (unit: nm) of the exponential function of
   the Yukawa potential. The Debye length of the homogeneous liquid will
   be used if the characteristic length is not specified. The dielectric
   constant for the Yukawa potential is defined in **-dielect-y**.

-  **-ccutoff**: the cutoff of closures, which is used in the closures
   like PLHNC.

-  **-rb** or **-Bohr-radius**: the minimal raiuds of an atom. By
   default it is 0.052911 (nm).

-  **-sd** or **-significant-digits**: the significant digits in the
   TS4S files, can be “float”, “double” or a number (significant digits
   in decimal)

-  **-closure-enhance-level** 1 or **-enhance-closure** 1: scale down
   the changes of SCF iterations with :math:`(1+h^2)^{\alpha/2}`, where
   :math:`\alpha` is the closure enhancement level. This option greatly
   helps the convergence of the self-consistent-field iterations and is
   turned on by default. Use **-no-closure-enhance** to turn this
   feature off.

-  **-bounded-to-ram**: don’t exceed the physical memory capacity when
   using memory. This feature is on by default, and use
   **-ignore-memory-capacity** to turn this feature off.

-  **-xvv-extend** 0: extend the **gvv** of solvent. If the input gvv
   contains the RDFs of :math:`2` nm, then -xvv-extend 5 will extend it
   to :math:`10` nm by filling the extended regions with :math:`1`. This
   option is helpful when you **gvv** is poor, but don’t expect too
   much.

.. |The screen output of running 3DRISMHI on alanine dipeptide| image:: rismhi3d-demo-screen.png
