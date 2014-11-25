__author__ = 'Vilius'

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


class Analyzer:
    def __init__(self, protein_path, search_output_file, blasta_file, mafft_output, program, database, entrez_query,
                 query_coverage):
        self.protein_path = protein_path
        self.record = self.load_protein()
        self.program = program
        self.database = database
        self.entrez_query = entrez_query
        self.queryCoverage = query_coverage
        self.search_output_file = search_output_file
        self.blasta_file = blasta_file
        self.mafft_output = mafft_output

    def load_protein(self):
        return SeqIO.read(open(self.protein_path), format="fasta")

    def blast(self):
        blast_records = self.get_blast_records()
        query_letters = blast_records.query_letters
        not_valid_alignments = []

        for alignment in blast_records.alignments:
            align_length = 0
            for hsp in alignment.hsps:
                align_length += hsp.align_length
            coverage = align_length / query_letters * 100
            if coverage < self.queryCoverage:
                not_valid_alignments.append(alignment)

        for alignment in not_valid_alignments:
            blast_records.alignments.remove(alignment)

        self.save_blast_record(blast_records)

    def get_blast_records(self):
        if not os.path.isfile(self.search_output_file):
            blast_output = NCBIWWW.qblast(self.program, self.database, self.record.seq,
                                          entrez_query=self.entrez_query, hitlist_size=500, expect=100.0)

            output_file = open(self.search_output_file, "w")
            output_file.write(blast_output.read())
            output_file.close()

        blast_file = open(self.search_output_file)
        return NCBIXML.read(blast_file)

    def save_blast_record(self, blast_records):
        fasta_list = []
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                fasta_list.append(">{0}\n{1}\n\n".format(alignment.title, hsp.sbjct))

        with open(self.blasta_file, "w") as tempFile:
            tempFile.writelines(fasta_list)

    def mafft(self, output="mafft_output.fasta"):
        os.system("mafft\mafft --localpair --quiet {0} > {1}".format(self.blasta_file, output))
        return

    def analyze(self, frag_len):
        all_seqs = []
        for r in SeqIO.parse(open(self.mafft_output), format="fasta"):
            all_seqs.append(r.seq)

        main_sequence = all_seqs[0]
        min_differences = 0
        min_index = 0
        max_differences = 0
        max_index = 0

        for i in range(0, len(main_sequence) - frag_len):
            differences = self.differences_in_fragment(main_sequence, all_seqs, i, frag_len)

            if differences > max_differences:
                min_differences = differences
                min_index = i
            if differences < min_differences:
                min_differences = differences
                max_index = i

        print("Mažiausiai panaši seka: {0}".format(main_sequence[min_index:(min_index + frag_len)]))
        print("Labiausiai panaši seka: {0}".format(main_sequence[max_index:(max_index + frag_len)]))
        return

    @staticmethod
    def differences_in_fragment(main_sequence, all_seqs, start_index, frag_len):
        differences = 0

        for i in range(1, len(all_seqs)):
            for j in range(start_index, start_index + frag_len):
                if all_seqs[i][j] != main_sequence[j]:
                    differences += 1

        return differences


def main():
    analyzer = Analyzer("data\\serum_albumin_preproprotein.fasta",
                        "data\\search_output.xml",
                        "data\\blast_output.fasta",
                        "data\\mafft_output.fasta"
                        "blastp",
                        "swissprot",
                        '"serum albumin"[Protein name] AND mammals[Organism]',
                        80)
    analyzer.blast()
    analyzer.mafft()
    analyzer.analyze(15)


main()