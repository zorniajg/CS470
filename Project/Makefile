default: dijkstra_serial par_dijkstra

dijkstra_serial: dijkstra_serial.c
	gcc -g -O2 --std=c99 -Wall -o dijkstra_serial dijkstra_serial.c

par_dijkstra: par_dijkstra.c
	mpicc -g -O2 --std=c99 -Wall -o par_dijkstra par_dijkstra.c

clean:
	rm -f dijkstra_serial par_dijkstra
