		
%.o : %.c  
	gcc -c $< -o $@ -I/sw/include/
	
main : main.o 
	gcc main.o -o main -ltiff -lm -L/sw/lib

default : main	
		
clean : 
	rm -f main.o main
