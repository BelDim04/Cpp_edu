#!/bin/bash

compile_debug(){
	echo -n -e "\e[36mSanitize check: \e[0m"
	cmd=$(make debug filein=$1 fileout=$2)
	res=$($cmd > /dev/null  2>&1)
	if [[ $res == "" ]]; then
		echo -e "\e[32mOK\e[0m"
	else
		echo -e "\n$res"
		exit -1
	fi
}

style_check(){
	echo -n -e "\e[36mCode-style check: \e[0m" 
	res=$(clang-format --dry-run --style=Google $1 > /dev/null 2>&1)
	if [[ $res == "" ]]; then
		echo -e "\e[32mOk\e[0m"
	else
		echo -e "\n$res"
		exit -1
	fi
}

compile_release(){
	$(make release filein=$1 fileout=$2 > /dev/null 2>&1)
}

compile_debug $@ && style_check $@ && compile_release $@ &&echo -e "\e[36mRun and Time: \e[0m" && time  ./$2
