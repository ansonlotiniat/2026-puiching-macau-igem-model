"""Read cached numeric cells from an XLSX worksheet using the standard library."""
import zipfile
import xml.etree.ElementTree as ET
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

def sheet_cells(path, sheet_name):
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ["".join(t.text or "" for t in node.findall(".//m:t", NS))
                      for node in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        relationships = {e.attrib["Id"]: e.attrib["Target"] for e in
                         ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))}
        sheets = ET.fromstring(archive.read("xl/workbook.xml")).find("m:sheets", NS)
        sheet = next(s for s in sheets if s.attrib["name"] == sheet_name)
        target = relationships[sheet.attrib[REL + "id"]]
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        result = {}
        for cell in ET.fromstring(archive.read(target)).findall(".//m:sheetData/m:row/m:c", NS):
            value = cell.find("m:v", NS)
            if value is None or value.text is None:
                continue
            if cell.attrib.get("t") == "s":
                result[cell.attrib["r"]] = shared[int(value.text)]
            elif cell.attrib.get("t") == "str":
                result[cell.attrib["r"]] = value.text
            else:
                result[cell.attrib["r"]] = float(value.text)
        return result
