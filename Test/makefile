#compilation pour les TP et pour lieer les libraries bien

# Options de compilation																													  v optimisation flag -Ofast desconseille de utiliser dans les codes paralele et les codes pas tres robust/teste
CFLAGS = -Wall -Wextra -Wfloat-equal -Wundef -Wshadow -Wcast-align  -Wcast-qual -Wunreachable-code -Wfloat-equal -Wunsafe-loop-optimizations #-O2 #-O 1 2 3 fast g s# -Wconversion  -fsanitize=address -Wstrict-prototypes#-fsanitize=thread -fsanitize=undefined -fstack-protector
#CFLAGS = -g																													# 			  ^  recomende de le laisser on -O2 car apres on voix pas de difference

# nom des librairies a lier     
LNAME = -lm 

# nom de l'executable
EXEC=EXE.exe


# https://riptutorial.com/makefile/example/21376/building-from-different-source-folders-to-different-target-folders

INCDIR  = ../Include
PROJDIR = $(CURDIR) # GNU make variable to the current directory
BUILDDIR= $(CURDIR)

VERBOSE = TRUE # to show or not show commands

#liste des fichiers sources
#SRC=  main.c lectureDonnee.c calculs.c
SRC= $(wildcard $(INCDIR)/*.c) $(wildcard *.c) #$(wildcard *.c) # $(wildcard */*.cpp) #si on peux regarde dans les dossiers 
#liste des fichiers objets
OBJ=  $(wildcard *.o) $(SRC:.c=.o)


# action declenchee par defaut
all: $(EXEC)

# compilation des fichiers sources en fichiers objets  
%.o: %.c
	$(CC) $(CFLAGS) -c $^

# compilation de l exectubale a partir des fichier sources
$(EXEC) : $(OBJ) 
	@echo why
	make wut
	chmod +x $(EXEC)
	@echo compilation terminee
# @echo "faites \"./cheerio.exe && make visualise\" pour visualiser le programme en temps reel et pour le sauvgarder en gif"
# les @ sont ici pour silencer les appels. Comme ca on voix que les resultats de ces fonctions mais pas lappel avel le echo 

# i need to put this otherwise it doesnt wait to compile files and says there is only one file 
wut:
	$(CC) -o  $(EXEC) $(wildcard *.o) $(LNAME)

test:
	echo $(SRC:.c=.o)
	echo $(wildcard *.o)

clean:
	rm -rf *.o
	ls



# #compilation pour les TP et pour lieer les libraries bien

# # Options de compilation																													  v optimisation flag -Ofast desconseille de utiliser dans les codes paralele et les codes pas tres robust/teste
# CFLAGS = -Wall -Wextra -Wfloat-equal -Wundef -Wshadow -Wcast-align  -Wcast-qual -Wunreachable-code -Wfloat-equal -Wunsafe-loop-optimizations #-O2 #-O 1 2 3 fast g s# -Wconversion  -fsanitize=address -Wstrict-prototypes#-fsanitize=thread -fsanitize=undefined -fstack-protector
# #CFLAGS = -g																													# 			  ^  recomende de le laisser on -O2 car apres on voix pas de difference

# # nom des librairies a lier     
# LNAME = -lm 

# # nom de l'executable
# EXEC=TP1.exe


# # https://riptutorial.com/makefile/example/21376/building-from-different-source-folders-to-different-target-folders

# INCDIR  = ../Include
# PROJDIR = $(CURDIR) # GNU make variable to the current directory
# BUILDDIR= $(CURDIR)

# VERBOSE = TRUE # to show or not show commands

# #liste des fichiers sources
# #SRC=  main.c lectureDonnee.c calculs.c
# SRC= $(wildcard $(INCDIR)/*.c) $(wildcard *.c) #$(wildcard *.c) # $(wildcard */*.cpp) #si on peux regarde dans les dossiers 
# #liste des fichiers objets
# OBJ=  $(wildcard *.o) $(SRC:.c=.o)


# # action declenchee par defaut
# all: $(EXEC)

# # compilation des fichiers sources en fichiers objets  
# %.o: %.c
# 	$(CC) $(CFLAGS) -c $^

# # compilation de l exectubale a partir des fichier sources
# $(EXEC) : $(OBJ)
# 	$(compile)
# 	chmod +x $(EXEC)
# 	@echo compilation terminee
# # @echo "faites \"./cheerio.exe && make visualise\" pour visualiser le programme en temps reel et pour le sauvgarder en gif"
# # les @ sont ici pour silencer les appels. Comme ca on voix que les resultats de ces fonctions mais pas lappel avel le echo 

# $(compile): $(wildcard *.o)
# 	$(CC)  -o $@ $^ $(LNAME)

# test:
# 	echo $(SRC:.c=.o)
# 	echo $(wildcard *.o)

# clean:
# 	rm -rf *.o
# 	ls

