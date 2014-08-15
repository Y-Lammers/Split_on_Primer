The Split_on_Primer.xml and Split_on_Primer_wrapper.py files have to be place in Galaxy's "tools" folder. 
Galaxy's "tool_conf.xml" file in the main folder needs to edited to include the HTS-barcode-checker, 
see: "http://wiki.galaxyproject.org/Admin/Tools/Adding_Tools" for details. In order for the tool to work 
the Split_on_Primer.xml (in src directory) needs to be added to the systems $PATH.
