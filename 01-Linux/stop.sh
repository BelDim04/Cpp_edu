#!/bin/bash

stop_server(){
	pid=$(cat server.pid 2>/dev/null)
	if [[ $pid == "" ]]; then
		echo no server
		exit -1
	fi
	kill $pid
	echo server stopped
}

del_pid(){
	rm -f server.pid
}

cd /home/judge
stop_server && del_pid
