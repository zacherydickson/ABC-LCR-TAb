#!/bin/bash

installDir=$1;
[ -z $installDir ] && installDir="$HOME/bin";
prog=$2;
[ -z $prog ] && prog="run_abc";

[ -d $installDir ] || { >&2 echo "Installation Dir ($installDir) doesn't exist"; exit 1; } ;

my_path="$(dirname -- ${BASH_SOURCE[0]})";

installFile="$installDir/$prog";
[ -s $installFile ] && { >&2 echo "$installFile Already exists"; exit 1; } ;

execPath=$(readlink -f ${my_path}/run_abc.sh);

echo -e "#!/bin/bash\n" >| $installFile;
echo -e "$execPath \$@" >> $installFile;
chmod +x $installFile;
