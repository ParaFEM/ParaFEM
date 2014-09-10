# Installing Extrae

1. Install dependencies: libxml2, binutils (bfd and libiberty), papi, zlib)
	See Section 3.5.2 BlueGene/Q in the User Guide: http://www.bsc.es/computer-sciences/performance-tools/trace-generation/extrae/extrae-user-guide
		Note: A static version of libxml2 is required.
		Note: The cross-compiler must be used for binutils.
2. Point the Linux-BGQ prebuilt version of Extrae to the dependencies.


LD_PRELOAD mechanism is not available on IBM machines so we must build with Extrae...



# Building w/ Extrae
1. Static linking is required (See joule_extrae.inc - line 67 for example)
2. Ensure all environment variables are set before building (See joule_extrae_vars.sh for example)




# Running w/ Extrae
1. Ensure you have a configuration XML file for Extrae (See extrae.xml for example)
2. Check the binary defined in the XML (See extrae.xml - line 66 for example)
3. Ensure the job script points to the XML file (See ll_trace.sh - line 3 and line 21)
4. Submit as usual.
