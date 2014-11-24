#Sis skrtipukas vykdo NCBI blast uzklausa ir atspausdina rastus atitikmenis bei  sekas
import Bio.Blast.NCBIWWW
import Bio.Blast.NCBIXML
enterez_query = '"serum albumin"[Protein name] AND mammals[Organism]'
ncbi = Bio.Blast.NCBIWWW.qblast(program="blastp" , database="swissprot", sequence="4502027", entrez_query=enterez_query,  hitlist_size=500,  expect=100.0) #seka, pagal kuria ieskoma  galima nurodyti ir kodu
# Nuskaitom XML formatu parsiustus duomenis. Patartina siuos duomenis issisaugoti i faila, ir pakartotinai kreiptis tik esant butinybei.
blast = Bio.Blast.NCBIXML.read(ncbi)
for sequence in blast.alignments:
    print('>%s'%sequence.title) # rasto atitikmens pavadinimas fasta formatu
    print(sequence.hsps[0].sbjct) # rasto atitikmens labiausiai patikimo sutampancio fragmento seka. Kiti maziau patikimi fragmentai [kuriu indeksai didesni nei 0] nedomina.
