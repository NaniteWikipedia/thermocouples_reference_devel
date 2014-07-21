
all: thermocouples_reference/source_OMEGA.py thermocouples_reference/source_NIST.py

thermocouples_reference/source_OMEGA.py: create_tables_OMEGA.py
	python $<  >  $@

thermocouples_reference/source_NIST.py: create_tables_NIST.py
	python $<  >  $@
