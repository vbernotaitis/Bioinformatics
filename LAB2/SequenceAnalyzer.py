__author__ = 'Vilimantas Bernotaitis'

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO


class SequenceAnalyzer:
    def __init__(self, protein_path, search_output_file, blasta_file, mafft_output, program, database, entrez_query, query_coverage):
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

        for alignment in blast_records.alignments:
            align_length = 0
            for hsp in alignment.hsps:
                align_length += hsp.align_length
            coverage = align_length / blast_records.query_letters * 100
            if coverage < self.queryCoverage:
                blast_records.alignments.remove(alignment)

        self.save_blast_record(blast_records)

    def get_blast_records(self):
        if not os.path.isfile(self.search_output_file):
            blast_output = NCBIWWW.qblast(self.program, self.database, self.record.seq, self.entrez_query, 500, 100.0)
            with open(self.search_output_file, "w") as tempFile:
                tempFile.write(blast_output.read())

        blast_file = open(self.search_output_file)
        return NCBIXML.read(blast_file)

    def save_blast_record(self, blast_records):
        fasta_list = []
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                fasta_list.append(">{0}\n{1}\n\n".format(alignment.title, hsp.sbjct))

        with open(self.blasta_file, "w") as tempFile:
           tempFile.writelines(fasta_list)

    def mafft(self):
        os.system("mafft\mafft --localpair --quiet {0} > {1}".format(self.blasta_file, self.mafft_output))
        return

    def analyze(self, frag_len):
        self.blast()
        self.mafft()

        all_seqs = []
        for r in SeqIO.parse(open(self.mafft_output), format="fasta"):
            all_seqs.append(r.seq)

        main_sequence = all_seqs[0]
        min_differences = len(main_sequence) * len(all_seqs)-1
        min_index = 0
        max_differences = 0
        max_index = 0

        for i in range(0, len(main_sequence) - frag_len):
            differences = self.differences_in_fragment(main_sequence, all_seqs, i, frag_len)

            if differences < min_differences:
                min_differences = differences
                min_index = i
            if differences > max_differences:
                max_differences = differences
                max_index = i

        print("Originaliausia seka: {0}".format(main_sequence[max_index:(max_index + frag_len)]))
        print("Index: {0}".format(max_index))
        print("Skirtumai: {0}".format(max_differences))
        for i in range(0, len(all_seqs)):
            print(all_seqs[i][max_index:max_index+frag_len])

        print("Labiausiai panaši seka: {0}".format(main_sequence[min_index:(min_index + frag_len)]))
        print("Index: {0}".format(min_index))
        print("Skirtumai: {0}".format(min_differences))
        for i in range(0, len(all_seqs)):
            print(all_seqs[i][min_index:min_index+frag_len])
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
    analyzer = SequenceAnalyzer("data\\serum_albumin_preproprotein.fasta",
                        "data\\search_output.xml",
                        "data\\blast_output.fasta",
                        "data\\mafft_output.fasta",
                        "blastp",
                        "swissprot",
                        '"serum albumin"[Protein name] AND mammals[Organism]',
                        80)

    analyzer.analyze(15)


main()