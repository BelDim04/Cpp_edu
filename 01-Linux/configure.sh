#!/bin/bash

create_user(){
	echo "Creating user (judge)..."
	res=$(useradd -m judge > /dev/null 2>&1)
	if [[ $res == "" ]]; then
		cd /home/judge
		echo User created
	else
		echo $res
		exit -1
	fi
}

download_artefact(){
	echo "Download artefact ($1)..."
	res=$(wget -O artefact.tar.gz $1 > /dev/null 2>&1)
	if [[ $res == "" ]]; then
		echo Artefact downloaded
	else 
		echo $res
		exit -1
	fi
}

tar_artefact(){
	echo "Tar artefact (artefact.tar.gz)..."
	res=$(tar -xf artefact.tar.gz > /dev/null 2>&1)
	if [[ $res == "" ]]; then
		echo OK
	else
		echo $res
		exit -1
	fi
}

files_judge_set(){
	chown -R judge /home/judge
	chmod u+x server
}

judge_path_update(){
	p=$(find /home/judge -name "server")
	P=$PATH:$p
	echo PATH=$P >> /home/judge/.bashrc
}

create_user && download_artefact $@ && tar_artefact && files_judge_set && judge_path_update
