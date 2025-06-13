# Except following variables, one should source thisroot.sh
unset LCIO
export LCIO=$PWD/../LCIO 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LCIO/lib

export DRUIDDIR=$HOME/opt/druid/Druid-2.4

export PATH=$PATH:$DRUIDDIR:$LCIO/bin:$DRUIDDIR/bin 

