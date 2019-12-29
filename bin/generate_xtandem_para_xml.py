import sys


def getParaItem(type, label, value):
	out = "<note type=\""+type+"\" label=\""+ label+"\""+">"+str(value)+"</note>\n"
	return(out)

def createTaxFile(db, out_file, db_name="tandem"):
	od = open(out_file,"w")
	od.write("<?xml version=\"1.0\"?>\n")
	od.write("<bioml label=\"x! taxon-to-file matching list\">\n")
	od.write(" <taxon label=\""+db_name+"\">\n")
	od.write("  <file format=\"peptide\" URL=\""+db+"\" />\n")
	od.write(" </taxon>\n")
	od.write("</bioml>\n")
	od.close()


def createInputXML(default_xml, tax_xml, ms_file, out_search_file, out_xml_file):
	of = open(out_xml_file,"w")

	of.write("<?xml version=\"1.0\"?>\n")
	of.write("<bioml>\n")

	of.write(getParaItem("input","list path, default parameters",default_xml))
	of.write(getParaItem("input","list path, taxonomy information",tax_xml))
	of.write(getParaItem("input","spectrum, path",ms_file))
	of.write(getParaItem("input","protein, taxon","tandem"))
	of.write(getParaItem("input","output, path",out_search_file))

	of.write("</bioml>\n")

	of.close()
	

## parameter for MS/MS searching
input_xml  = sys.argv[1]
## MS/MS data file path
ms_path    = sys.argv[2]
## Database file path
db_path    = sys.argv[3]
## Output file prefix
out_prefix = sys.argv[4]
## Output xml file
out_xml    = sys.argv[5]


tax_file = str(out_prefix) + "_tax.xml"
createTaxFile(db_path, tax_file)
out_search_file = str(out_prefix) + ".xml"
createInputXML(input_xml, tax_file, ms_path, out_search_file, out_xml)





