#!/bin/bash
# usage: ./run_eik.sh MESH_um NAME X_um Y_um Z_um OTHER_ARGS...   (MESH_um without extension)
# The eikonal domain (phys_region) and the gregions with the velocities are passed by
# run_sims.py. The FIM stops after 1000 iterations: the LATs stop changing after
# about 120, even in the runs that do not go below fim.tol.
MESH=$1; NAME=$2; X=$3; Y=$4; Z=$5; shift 5
openCARP \
 -meshname $MESH -simID run_$NAME \
 -dream.solve 0 -dream.fim.tol 0.001 -dream.fim.max_iter 1000 \
 -dream.output.atfile act -tend 500 \
 -num_stim 1 -stim[0].elec.geom_type 1 \
 -stim[0].elec.p0[0] $X -stim[0].elec.p0[1] $Y -stim[0].elec.p0[2] $Z \
 -stim[0].elec.radius 3000 -stim[0].ptcl.start 0.0 \
 "$@"
