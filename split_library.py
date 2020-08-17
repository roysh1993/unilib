from Bio import SeqIO, Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import pandas as pd
import os
import sys
import time
lexA_sec_barcode_count = 0
universal_sec_barcode_count = 0
lexA_first_barcode_count = 0
universal_first_barcode_count = 0
total_reads_count = 0
junk_counter = 0


def stupid_split_library(fastq_file,barcode_df,universal_out,lexA_out):
    global lexA_sec_barcode_count,universal_sec_barcode_count,lexA_first_barcode_count,\
        universal_first_barcode_count,total_reads_count

    univeral_first_barcodes = list(barcode_df["first_barcode_universal"].dropna())
    lexA_first_barcodes = list(barcode_df["first_barcode_lexA "].dropna())
    lexA_buffer = []
    universal_buffer = []
    start = time.time()
    with open(fastq_file, "rU") as handle, open(universal_out, "w") as universal_handle, \
            open(lexA_out, "w") as lexA_handle:
        for record in SeqIO.parse(handle, "fastq"):

            single_read = str(record.format("fastq")).split("\n")
            dna_string = single_read[1]
            first_barcode = dna_string[:7]

            # first barcode is of universal library
            if first_barcode in univeral_first_barcodes:
                universal_buffer.append(record)
                universal_first_barcode_count += 1

            if first_barcode in lexA_first_barcodes:
                lexA_buffer.append(record)
                lexA_first_barcode_count += 1

            if len(lexA_buffer) == 100000:
                SeqIO.write(lexA_buffer, lexA_handle, "fastq")
                lexA_buffer = []

            if len(universal_buffer) == 100000:
                SeqIO.write(universal_buffer, universal_handle, "fastq")
                universal_buffer = []

            if total_reads_count % 1000000 == 0:
                print("after ",total_reads_count," iterations, total time:",time.time()- start)

            total_reads_count += 1
        SeqIO.write(lexA_buffer, lexA_handle, "fastq")
        SeqIO.write(universal_buffer, universal_handle, "fastq")
        print("total time:", time.time() - start)


def efficient_split_library(fastq_file,barcode_df,universal_out,lexA_out,junk_out):
    global lexA_sec_barcode_count,universal_sec_barcode_count,lexA_first_barcode_count,\
        universal_first_barcode_count,total_reads_count, junk_counter
    print(barcode_df["sec_barcode_universal"].isna().sum())
    print(barcode_df["sec_barcode_lexA"].isna().sum())
    univeral_first_barcodes = barcode_df["first_barcode_universal"].dropna().values
    lexA_first_barcodes = barcode_df["first_barcode_lexA "].dropna().values
    lexA_buffer = []
    junk_buffer = []
    universal_buffer = []

#     create dictionery
    universal_hash = set(barcode_df["sec_barcode_universal"].dropna())
    lexA_hash = set(barcode_df["sec_barcode_lexA"].dropna())

    start = time.time()
    with open(fastq_file, "rU") as handle, open(universal_out, "w") as universal_handle, \
            open(lexA_out, "w") as lexA_handle,open(junk_out, "w") as junk_handle:
        for record in SeqIO.parse(handle, "fastq"):

            single_read = str(record.format("fastq")).split("\n")
            dna_string = single_read[1]
            first_barcode = dna_string[:7]

            # first barcode is of universal library
            if first_barcode in univeral_first_barcodes:
                universal_first_barcode_count += 1
                #     test second barcode
                second_barcode = dna_string[39:54]
                if second_barcode in universal_hash:
                    universal_sec_barcode_count += 1
                    universal_buffer.append(record)

            elif first_barcode in lexA_first_barcodes:
                lexA_first_barcode_count += 1
                # test second barcode
                second_barcode = dna_string[41:56]
                #     test second barcode
                if second_barcode in lexA_hash:
                    lexA_sec_barcode_count += 1
                    lexA_buffer.append(record)
            else:
                junk_buffer.append(record)
                junk_counter += 1




            total_reads_count += 1
#                 if buffer full append
            if len(lexA_buffer) == 100000:
                SeqIO.write(lexA_buffer, lexA_handle, "fastq")
                lexA_buffer = []
            if len(universal_buffer) == 100000:
                SeqIO.write(universal_buffer,universal_handle,"fastq")
                universal_buffer = []
            if len(junk_buffer) == 100000:
                SeqIO.write(junk_buffer,junk_handle,"fastq")
                junk_buffer = []

            if total_reads_count % 1000000 == 0:
                print("after ",total_reads_count," iterations, total time:",(time.time() - start)/360)

        # write rest of buffer
        SeqIO.write(lexA_buffer, lexA_handle, "fastq")
        SeqIO.write(universal_buffer, universal_handle, "fastq")
        SeqIO.write(junk_buffer, junk_handle, "fastq")

        # save stats to file
        stats_file = os.path.splitext(fastq_file)[0] + "split_stats.txt"
        with open(stats_file,'w') as f:
            print(f'Time: {start -time.time()}')
            print(f'total number of reads: {total_reads_count}')
            f.write(f'total number of reads: {total_reads_count}\n')
            print(f'lexA first barcode count: {lexA_first_barcode_count}')
            f.write(f'lexA first barcode count: {lexA_first_barcode_count}\n')
            print(f'lexA second barcode count: {lexA_sec_barcode_count}')
            f.write(f'lexA second barcode count: {lexA_sec_barcode_count}\n')
            print(f'universal first barcode count: {universal_first_barcode_count}')
            f.write(f'universal first barcode count: {universal_first_barcode_count}\n')
            print(f'universal second barcode count: {universal_sec_barcode_count}')
            f.write(f'universal second barcode count: {universal_sec_barcode_count}\n')
            print(f'junk count: {junk_counter}')
            f.write(f'junk count: {junk_counter}')
        print("DONE")


def split_library(fastq_file,barcode_df,universal_out,lexA_out):
    global lexA_sec_barcode_count,universal_sec_barcode_count,lexA_first_barcode_count,\
        universal_first_barcode_count,total_reads_count

    univeral_first_barcodes = barcode_df["first_barcode_universal"].values
    lexA_first_barcodes = barcode_df["first_barcode_lexA "].values

    universal_sec_barcodes = barcode_df["sec_barcode_universal"].values
    lexA_sec_barcodes = barcode_df["sec_barcode_lexA"].values
    start = time.time()
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

            if total_reads_count % 1000 == 0:
                print("after ",total_reads_count," iterations, total time:",(time.time() - start)/360)

            total_reads_count += 1


            # check first barcode return
    print(f'Time: {start - time.time()}')
    print(f'total number of reads: {total_reads_count}')
    print(f'lexA first barcode count: {lexA_first_barcode_count}')
    print(f'lexA second barcode count: {lexA_sec_barcode_count}')
    print(f'universal first barcode count: {universal_first_barcode_count}')
    print(f'universal second barcode count: {universal_sec_barcode_count}')







if __name__ == '__main__':

    barcode_df = pd.read_excel("merged_barcodes.xlsx")
    efficient_split_library(sys.argv[1],barcode_df,sys.argv[2],sys.argv[3],sys.argv[4])




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
