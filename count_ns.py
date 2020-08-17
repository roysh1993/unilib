from Bio import SeqIO, Seq
import os
import sys

def count_ns_in_barcode(file_name,barcode_start,barcode_end):
    n_0 = 0
    n_1 = 0
    n_2 = 0
    n_more = 0
    total_count = 0
    output_file = os.path.splitext(file_name)[0] +"n_count.txt"

    with open(file_name, "rU") as handle:

        for record in SeqIO.parse(handle, "fastq"):
            single_read = str(record.format("fastq")).split("\n")
            dna_string = single_read[1]
            first_barcode = dna_string[:7]
            number_of_n = first_barcode.count("N")

            #switch case
            if number_of_n == 0:
                n_0 +=1
            elif number_of_n == 1:
                n_1 += 1
            elif number_of_n == 2:
                n_2 += 1
            else:
                n_more += 1

            total_count += 1

        #write results to file

        with open(output_file,'w') as out:
            out.write("N count of file {} \n".format(file_name))
            out.write("Total number of reads: {} \n".format(total_count))
            out.write("Number of Barcodes with 0 N's: {} \n".format(n_0))
            out.write("Number of Barcodes with 1 N's: {} \n".format(n_1))
            out.write("Number of Barcodes with 2 N's: {} \n".format(n_2))
            out.write("Number of Barcodes with more then 2 N's: {} \n".format(n_more))

        print("done")
        print("result output on file: ",output_file)


if __name__ == '__main__':
    count_ns_in_barcode(sys.argv, 0, 7)



