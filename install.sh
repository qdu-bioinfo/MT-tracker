#!/bin/bash
# MT-tracker Installer
# Bioinformatics Group, Qingdao University
# Updated at Jan. 16, 2025
# Updated by Wenjie Zhu, Xiaoquan Su

#########################################
# 1. Determine the shell configuration file
#########################################
if [[ $SHELL = '/bin/zsh' ]]; then
    CONFIG_FILE=~/.zshrc
    if [ ! -f "$CONFIG_FILE" ]; then
        CONFIG_FILE=~/.zsh_profile
        if [ ! -f "$CONFIG_FILE" ]; then
            touch $CONFIG_FILE
        fi
    fi
else
    CONFIG_FILE=~/.bashrc
    if [ ! -f "$CONFIG_FILE" ]; then
        CONFIG_FILE=~/.bash_profile
        if [ ! -f "$CONFIG_FILE" ]; then
            touch $CONFIG_FILE
        fi
    fi
fi

#########################################
# 2. Set installation path and system version
#########################################
MT_PATH=`pwd`
Sys_ver=`uname`

#########################################
# 3. Check for existing MTTRACKER environment variable and PATH settings in the config file
#########################################
Old_MT=`grep "export MTTRACKER" $CONFIG_FILE | awk -F '=' '{print $1}'`
Old_PATH=`grep "MTTRACKER/bin" $CONFIG_FILE | sed 's/\(.\).*/\1/' | awk '{if($1!="#"){print "True";}}'`
Disable_Prefix="####DisabledbyMTtracker####"

echo "**MT-tracker Installation**"
echo "**Version 1.0**"
echo "**Updated at Jan. 16, 2025**"

#########################################
# 4. Build the source code if a Makefile exists
#########################################
if [ -f "Makefile" ]; then
    echo -e "\n**MT-tracker src package build**"
    make
    echo -e "\n**Build Complete**"
else
    echo -e "\n**MT-tracker bin package**"
fi

#########################################
# 5. Configure the MTTRACKER environment variable in the config file
#########################################
if [ "$Old_MT" != "" ]; then
    Current_MT=`grep ^export\ MTTRACKER $CONFIG_FILE | awk -F '=' '{print $2}'`
    if [ "$Current_MT" != "$MT_PATH" ]; then
        if [ "$Sys_ver" = "Darwin" ]; then
            sed -i "" "s/^export\ MTTRACKER/$Disable_Prefix\ &/g" $CONFIG_FILE
            sed -i "" -e "`grep -n "$Disable_Prefix" $CONFIG_FILE | cut -d ":" -f 1 | head -1` a\\
export MTTRACKER=$MT_PATH
" $CONFIG_FILE
        else
            sed -i "s/^export\ MTTRACKER/$Disable_Prefix\ &/g" $CONFIG_FILE
            sed -i "/$Disable_Prefix\ export\ MTTRACKER/a export MTTRACKER=$MT_PATH" $CONFIG_FILE
        fi
    fi
else
    echo "export MTTRACKER=$MT_PATH" >> $CONFIG_FILE
fi

#########################################
# 6. Append MTTRACKER/bin to PATH if not already present
#########################################
if [ "$Old_PATH" = "" ]; then
    echo 'export PATH="$PATH:$MTTRACKER/bin"' >> $CONFIG_FILE
fi

#########################################
# 7. Source the configuration file to load the new environment variables
#########################################
source $CONFIG_FILE
echo -e "\n**Environment Variables Configuration Complete**"

#########################################
# 8. Final installation message
#########################################
echo -e "\n**MT-tracker Installation Complete**"
echo -e "\n**An example dataset with demo script is available in \"example\"**"
