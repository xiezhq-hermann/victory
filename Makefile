include make_pgi.inc

DIR := ${PWD}

ctqmc: 

	@sed -i '3s|.*|DIR         = '${DIR}'|' make_pgi.inc
	@mkdir bin

	@echo '---------------- DFFT ----------------------'
	(cd $(DIR)/lib/dfftpack; make; mv libdfftpack.a ../libdfftpack.a)
	@echo '---------------- SRC_MOD ----------------------'
	(cd $(DIR)/src/SRC_MOD; make)
	@echo
	@echo '---------------- SRC_PA --------------------'
	(cd $(DIR)/src/SRC_PA; make)
	@echo '---------------- MOVE EXECUTABLE FILE -----------------'
	(mv $(DIR)/bin/victory $(DIR)/out)


clean:
	(cd $(DIR)/src/SRC_MOD; make clean)
	(cd $(DIR)/src/SRC_PA; make clean)
	(cd $(DIR)/lib/dfftpack; make clean)
	rm -f lib/lib.a
	rm -f lib/libdfftpack.a
	rm -r bin

