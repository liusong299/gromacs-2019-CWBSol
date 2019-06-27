/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cmath>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include <unistd.h>
#include "../3drismhi/rismhi3d.cpp"

//int main_3drism(int argc, char *argv[]);

int gmx_3drism(int argc, char *argv[])
{
    const char *desc[] = {
        "[THISMODULE] The 3D reference interaction site model (3DRISM) is a powerful tool to study the thermodynamic and structural properties of liquids.",
        "However, for hydrophobic solutes,",
        "the inhomogeneity of the solvent density around them poses a great challenge to the 3DRISM theory.",
        "To address this issue, we have previously introduced the hydrophobic-induced density inhomogeneity theory (HI) for purely hydrophobic solutes.",
        "To further consider the complex hydrophobic solutes containing partial charges,",
        "here we propose the D2MSA closure to incorporate the short-range and long-range interactions with the D2 closure and the mean spherical approximation, respectively.",
        "We demonstrate that our new theory can compute the solvent distributions around real hydrophobic solutes in water and complex organic solvents that agree well with the explicit solvent molecular dynamics simulations."};

    if (argc > 1)
    {
        //argv[0] = "/Users/stephen/projects_local/gromacs-2019/src/gromacs/rismhi3d/rismhi3d";
        //execv(argv[0], argv);

    }

    return main_3drism(argc, argv);
}


