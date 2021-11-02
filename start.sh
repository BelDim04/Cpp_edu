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
start_server && write_pid && healthcheck &
