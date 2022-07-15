debug:
	g++ -g -DDEBUG -fsanitize=address -fsanitize=leak -fsanitize=undefined -o $(fileout) $(filein)
release:
	g++ -O3 -o $(fileout) $(filein)
