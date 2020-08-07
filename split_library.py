from Bio import SeqIO, Seq
import numpy as np
import pandas as pd
import time
lexA_sec_barcode_count = 0
universal_sec_barcode_count = 0
lexA_first_barcode_count = 0
universal_first_barcode_count = 0
total_reads_count = 0




def split_library(fastq_file,barcode_df,universal_out,lexA_out):
    global lexA_sec_barcode_count,universal_sec_barcode_count,lexA_first_barcode_count,\
        universal_first_barcode_count,total_reads_count

    univeral_first_barcodes = barcode_df["first_barcode_universal"].values
    lexA_first_barcodes = barcode_df["first_barcode_lexA "].values

    universal_sec_barcodes = barcode_df["sec_barcode_universal"].values
    lexA_sec_barcodes = barcode_df["sec_barcode_lexA "].values
    with open(fastq_file, "rU") as handle, open(universal_out, "w") as universal_handle,\
            open(lexA_out, "w") as lexA_handle:

        for record in SeqIO.parse(handle, "fastq"):
            single_read = str(record.format("fastq")).split("\n")
            dna_string = single_read[1]
            first_barcode = dna_string[:7]
            # first barcode is of universal library
            if first_barcode in univeral_first_barcodes:
                universal_first_barcode_count += 1
            #     test second barcode
                second_barcode = dna_string[39:54]
                if second_barcode in universal_sec_barcodes:
                    universal_sec_barcode_count += 1
            #         write to universal fastq
                    SeqIO.write(record,universal_handle , "fastq")

            if first_barcode in lexA_first_barcodes:
                lexA_first_barcode_count += 1
                #test second barcode
                second_barcode = dna_string[41:56]
            #     test second barcode
                if second_barcode in lexA_sec_barcodes:
                    lexA_sec_barcode_count += 1
                    SeqIO.write(record, lexA_handle, "fastq")

            total_reads_count += 1

            # check first barcode    return
    print(f'total number of reads: {total_reads_count}')
    print(f'lexA first barcode count: {lexA_first_barcode_count}')
    print(f'lexA second barcode count: {lexA_sec_barcode_count}')
    print(f'universal first barcode count: {universal_first_barcode_count}')
    print(f'universal second barcode count: {universal_sec_barcode_count}')

    return





if __name__ == '__main__':
    bardoce_df = pd.read_excel("merged_barcodes.xlsx")

    split_library("random_samp.fastq",bardoce_df,"uni_test.fastq","lexA_test.fastq")
    # a=5
    # print(bardoce_fastq.columns)
    # universal_barcodes = bardoce_fastq["sec_barcode_universal"].values
    # first_barcode = bardoce_fastq["first_barcode_lexA "].values
    #
    # s = time.time()
    # true_val = "GAAGGTTTCGGAGCC" in universal_barcodes
    # s = time.time() - s
    # print(true_val,"time: ",s)
    # s = time.time()
    # true_val = "GCATCTC" in first_barcode
    # s = time.time() - s
    # print(true_val, "time: ", s)
