# use this Makefile to submit in parallel. Adjust it and run:
# make -j10 SAMPLES=p1328/simplifiedModel_wA_slep.txt

# defines the tag of this processing
MYTAG=n0123

# set the username of the grid samples, by default it is your username on this machine
GRIDUSER=$(USER)

list = $(shell awk 'BEGIN{FS="/";ORS=" "}{print $$1}' $(SAMPLES))

all: $(list:%=%_$(MYTAG).prun)

%.prun:
	./submit.py mc -t $(MYTAG) -f $(SAMPLES) -p `basename $@ _$(MYTAG).prun` --destSE=TRIUMF-LCG2_LOCALGROUPDISK --nickname $(GRIDUSER) &> $@


