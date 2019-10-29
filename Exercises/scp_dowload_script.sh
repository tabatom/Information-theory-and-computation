scp tabarelt@spiro.fisica.unipd.it:/home/tabarelt/Quantum_Information/$1/$2/* $2
# Deve essere gestito meglio: voglio rinominare appositamente tutti gli eseguibili perché il compilatore fortran sulla spiro è vecchio rispetto a quello locale
# mv $1/$2/Ex*.x $1/$2/Ex*_old.x

# Some hints to use grep to get all executable files (named *.x)
# ls -lR | grep "\.x"
