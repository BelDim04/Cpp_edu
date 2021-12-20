#!/bin/bash

source healthcheck.sh

start_server(){
	./server &
	echo server started...
}

write_pid(){
	pid=$(pidof server)
	touch server.pid
	> server.pid
	echo $pid > server.pid
}

cd /home/judge
pid=$(cat server.pid  2>/dev/null)
if [[ $pid == "" ]]; then
	start_server && write_pid && healthcheck &
else
	echo server already exist
fi
