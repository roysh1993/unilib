from Bio import SeqIO, Seq
import sys
import json
# from scipy.io import loadmat

forward_primer = "TTAGAAATATCCGCAGCGCG"
reversed_primer = "CGCGCTGCGGATATTTCTAA"


g_f_counter = 0
g_r_counter = 0


def check_forward_primer(read):
    global g_f_counter
    result =  read[7:27] == forward_primer
    if result:
        g_f_counter += 1
    return result


def check_reverse_primer(read):
    global g_r_counter
    length = len(read)
    a = read[length-27:length -7]
    c = len(a)
    b = a[::-1]

    result = read[length-27:length -7]== reversed_primer
    if result:
        g_r_counter += 1

    return result


def covert_to_forward(fastq_input,fastq_output):

    with open(fastq_input, "rU") as handle, open(fastq_output, "w") as output_handle:

        for record in SeqIO.parse(handle, "fastq"):
            single_read = str(record.format("fastq")).split("\n")
            description = single_read[0]
            dna_string = single_read[1]
            quality_score = single_read[3]

            if check_forward_primer(dna_string):
                SeqIO.write(record,output_handle,"fastq")
            elif check_reverse_primer(dna_string):
                record = record.reverse_complement(id=True,description=True,name=True,)
                SeqIO.write(record,output_handle,"fastq")


# def load_mat_to_json(filename):
#     f_name = "data/UNILIB_post_E_coli_plasmid/ListStrcut60000.mat"
#     data = loadmat(f_name)
#     # j_s = json.dumps(data)
#     data = data["List"]
#     data2 = data[0][0][0]
#
#     a= 5
#
# def try_func():
#     f_name = "data/UNILIB_post_E_coli_plasmid/ListStrcut60000.mat"
#     matdata = loadmat(f_name)
#     # First, we'll combine the data sets for each unique combination of instruments, station, and depth.
#     matdata2 = {}  # Creates an empty dictionary
#
#     for k, v in matdata.iteritems():  # Here's way to loop over key,value pairs of dictionaries.
#         if k.startswith('__'):  # Skip the special keys
#             continue  # Continue goes to the next loop iteration.
#         param, station, date, depth = k.split(
#             '_')  # This splits up the string, assigning each of the parts to different variables.
#
#         # The if statements here just check to make sure that level of the dictionary has already been created.
#         # If all levels are present, we can just assign a new key.  If not, we have to create a new dictionary
#         # at that level.
#         if station in matdata2:
#             if param in matdata2[station]:
#                 if depth in matdata2[station][param]:
#                     matdata2[station][param][depth][date] = v  # We don't need to check depth b/c we're assigning it.
#                 else:
#                     matdata2[station][param][depth] = {date: v}
#             else:
#                 matdata2[station][param] = {depth: {date: v}}
#         else:
#             matdata2[station] = {param: {depth: {date: v}}}

if __name__ == '__main__':
    covert_to_forward("temp.fastq","blabla.fastq")
    print(g_r_counter)
    print(g_f_counter)

    # print(len("TTAGAAATATCCGCAGCGCG"))