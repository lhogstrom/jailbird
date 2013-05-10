#!/bin/bash

WORKINGDIR=/xchip/cogs/stacks/deprecated/STK010_T2D
OUTDIR=/xchip/cogs/hogstrom/analysis/T2D/gct

mkdir $OUTDIR
/xchip/cogs/cflynn/python_tools/bin/gct_13_to_12 -i $WORKINGDIR/STK010_QNORM_n756x978.gct -o $OUTDIR -f STK010_QNORM_n756x978_v12.gct

