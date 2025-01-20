###MT-tracker installer
###Bioinformatics Group, Qingdao University
###Updated at Jan. 16, 2025
###Updated by Wenjie Zhu,  Xiaoquan Su
#!/bin/bash

###Check Parallel-Meta environment variable###
if [ $ParallelMETA ]
    then
    echo -e "\n**Parallel-Meta is already installed**"
else
    echo -e "\n**Please install latest version of Parallel-Meta**"
    exit
fi

###Build source code for src package###
echo -e "\n**MT-tracker Installation**"
echo -e "\n**MT-tracker src package build**"
make
echo -e "\n**Build Complete**"

###Plugin installation###
echo -e "\n**Plugin Installation**"
cp bin/PM-Mt-tracker $ParallelMETA/bin
tar -xzvf metaphlan2.tar.gz
cp -rf mip_16s $ParallelMETA/databases

##Check database configuration##
Check_db_config=`grep -c "MetaPhlAn2" $ParallelMETA/databases/db.config`
if [ "$Check_db_config" = "0" ]
    then
    cp $ParallelMETA/databases/db.config $ParallelMETA/databases/.db.config.bk
    cat $ParallelMETA/databases/.db.config.bk db.config > $ParallelMETA/databases/db.config
fi

echo -e "\n**MT-tracker Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"