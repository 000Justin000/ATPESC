#! /bin/bash
echo ""
echo "4096 x 4096"
runjob --np 1024 -p 16 --block $COBALT_PARTNAME : mlife2d -x 4096 -y 4096 -i 100
