#!/bin/bash

check(){
	pid=$(cat server.pid  2>/dev/null)
	if [[ $pid == "" ]]; then
		echo "Server stopped"
		exit 0
	fi
	s=$(ps $pid)
	echo $s
}

healthcheck(){
	cd /home/judge
	touch server.log
	for ((;;))
	do
		r=$(check)
		echo $r >> server.log
		if [[ $r == "Server stopped" ]]; then
			break 1
		fi
		sleep 5
	done
}


