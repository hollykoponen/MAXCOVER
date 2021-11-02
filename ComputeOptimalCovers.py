import os
import subprocess
import gzip
import shutil
from random import sample
from Bio import SeqIO

# -- Globals -- #
gz_dir = "./test_data/GZ/"
faa_dir = "./test_data/FAA/"
data_dir = "./test_data/data_files/"
test_result_dir = "./test_data/test_results/"


# -- Functions -- #
def make_dirs(path):
    if not os.path.exists(path):
        os.makedirs(path)

# Unzips files of format {name}.faa.gz from gz_dir into faa_dir
def unzip_gz(filename):
    with gzip.open(os.path.join(gz_dir, filename), 'rb') as gz_in:
        out_name = os.path.join(faa_dir, os.path.splitext(filename)[0])
        with open(out_name, 'wb') as gz_out:
            shutil.copyfileobj(gz_in, gz_out)

# Sets Fasta_ID value
def compute_fasta_ID(fasta):
    if "|" in fasta.description:
        fasta_ID = fasta.description.split("|")[1]
    else:
        fasta_ID = fasta.description.split(" ")[0]
    return fasta_ID

# Clean data from min-max length, no X, and unique sequences. Superfolder is original Fasta filename
def clean_data(filename):
    fasta_sequences = SeqIO.parse(open(os.path.join(faa_dir,filename)), 'fasta')
    ID_Set = set()
    for fasta in fasta_sequences:
        fasta_ID = compute_fasta_ID(fasta)
        if ((len(fasta.seq) >= 20) 
                and (len(fasta.seq) <= 100000) 
                and (str(fasta.seq).find('X') == -1) 
                and (fasta_ID not in ID_Set)):
            ID_Set.add(fasta_ID)
    return ID_Set

# Based on sample sequence id list, extract these sequences from the fasta file
def extract_sample_seq(filename, sample_dir, SampleID_Set, include_header):
    fasta_sequences = SeqIO.parse(open(os.path.join(faa_dir,filename)), 'fasta')
    ctr = 0
    for fasta in fasta_sequences:
        fasta_ID = compute_fasta_ID(fasta)
        if fasta_ID in SampleID_Set:
            ctr += 1
            out_path = os.path.join(sample_dir, fasta_ID)
            with open("{}_{}.txt".format(out_path, ctr), 'w') as output_f:
                if include_header in ("Y", "y"): 
                    output_f.write(str("> {}\n".format(fasta.description)))            
                output_f.write(str(fasta.seq))

# Call the optimal cover C++ script
def Compute_OC(source, target, compute_OC, compute_top_ten, compute_repeatmatches, repeat_size):
    result_directory = os.path.join(test_result_dir, target)
    make_dirs(result_directory)
    print("\n{}".format(target))
    for file in os.listdir(source):
        result_file = os.path.join(result_directory, 'test_results_{}'.format(file))
        
        if file.endswith(".txt"):
            with open(result_file, "w") as f:
                abs_path = os.path.abspath(os.path.join(source, file))
                
                cmd = ["optimal_cover_script/bin/Release/optimal_covers", 
                    str(compute_OC),
                    str(compute_top_ten),
                    str(compute_repeatmatches),
                    str(repeat_size), 
                    abs_path]

                sp = subprocess.run(cmd, capture_output=True)     
                
                                
                f.writelines(sp.stdout.decode('utf-8'))
                fasta_sequences = SeqIO.parse(open(os.path.join(faa_dir,target+".faa")), 'fasta')
                for fasta in fasta_sequences:
                    fasta_ID = compute_fasta_ID(fasta)
                    if file.split("_")[0] == fasta_ID:
                        f.writelines(fasta.description)
                        
        print("{}".format(file))






# -- Main -- #
print("Script start")
make_dirs(gz_dir)
make_dirs(faa_dir)
make_dirs(data_dir)
make_dirs(test_result_dir)
Unzip = input("Unzip .faa.gz files? (Y/N) ")



# Unzip
while not os.listdir(gz_dir) and (Unzip in ["Y", "y"]):
    Unzip = input("No zip files in ./test_data/GZ/ folder. Add .faa.gzip file to this folder.\nUnzip .faa.gz files? (Y/N) ")

if Unzip in ("Y", "y"):
    for filename in sorted(os.listdir(gz_dir)):   
        if filename.endswith(".faa.gz"):
            unzip_gz(filename)
            print("Unzipped {}".format(filename)) 



# Split the Fasta File (all sequences or sample)
ExtractMaxData = input("Extract maximum number of clean data from dataset? (Y/N) ")

if ExtractMaxData in ("N", "n"):
    SampleSeq = input("Sample sequences? (Y/N) ")

if (ExtractMaxData in ("Y", "y")) or (SampleSeq in ("Y", "y")):
    include_header = input("Include the file header? i.e. fasta sequence description (Y/N) ")
    for filename in os.listdir(faa_dir):
        if filename.endswith(".faa"):
            sample_dir = os.path.join(data_dir, os.path.splitext(filename)[0])
            if os.path.exists(sample_dir):
                shutil.rmtree(sample_dir)
            make_dirs(sample_dir)
            ID_set = list(clean_data(filename))
            if ExtractMaxData in ("Y", "y"):
                SampleSize = len(ID_set) 
            elif SampleSeq in ("Y", "y"):
                SampleSize = int(input("Sample Size? (int) ")) 
            while int(len(ID_set)/SampleSize) < 1:
                SampleSize = int(input("Sample size too large. Max value for this dataset is {}.\nNew Sample Size? (int) ".format(len(ID_set))))
            SampleID_Set = sample(ID_set, SampleSize)
            extract_sample_seq(filename, sample_dir, SampleID_Set, include_header)
            print("Sampled {} sequences from {}".format(SampleSize, filename))



# Determine if user would like to compute all the optimal covers in the data files
WantToComputeOC = input("Compute Optimal Covers for all folders in /data_files/? (Y/N) ")

if WantToComputeOC in ["Y", "y"]:
    compute_OC = 1 #True
    compute_top_ten = 0 #False 
    compute_repeatmatches = 0 #False
    repeat_size = 1 #Default value - not used b/c only required for NE computation
    print("\nComputed Optimal Covers for:")
    for subfolder_name in sorted(os.listdir(data_dir)):
        source = os.path.join(data_dir, subfolder_name)
        Compute_OC(source, subfolder_name, compute_OC, compute_top_ten, compute_repeatmatches, repeat_size)

elif WantToComputeOC in ["N", "n"]:
    # Get Optimal Cover from a single folder?
    ComputeOC_SingleFolder = input("Do you want to compute Optimal Covers for a single folder in /data_files/? (Y/N)" )
    if ComputeOC_SingleFolder in ["Y", "y"]:
        compute_OC = 1 #True
        compute_top_ten = 0 #False 
        compute_repeatmatches = 0 #False
        repeat_size = 1 #Default value - not used b/c only required for NE computation
        OnlyOneFolder = input("Subfolder name? ")
        while not os.path.exists(os.path.join(data_dir, OnlyOneFolder)) and (OnlyOneFolder != "exit"):
             OnlyOneFolder = input("Directory does not exist.\n{}\nSubfolder name? (or type 'exit') ".format(os.path.join(data_dir, OnlyOneFolder)))
        if OnlyOneFolder == "exit":
            exit()
        source = os.path.join(data_dir, OnlyOneFolder)
        print("\nComputed Optimal Covers for:")
        Compute_OC(source, OnlyOneFolder, compute_OC, compute_top_ten, compute_repeatmatches, repeat_size)
     #TODO: Add Compute_top_ten into path
     #TODO: Put all Python Code in to C/C++ and rename MAXCOVER


# Determine if user would like to compute all repeat matches
Compute_RepeatMatches = input("Compute Repeat Matches? (Y/N) ")

if Compute_RepeatMatches in ("Y", "y"):
    repeat_size = int(input("How long should the substring length be for the repeat? (int, 1+) "))
    while (type(repeat_size) is not int) or (not repeat_size >= 1): 
        repeat_size = input("This is not an integer number greater than or equal to 1. Please enter an integer value for the substring length: ")
    compute_OC = 0 #False
    compute_repeatmatches = 1 #True
    print("\nComputed Repeat Matches for:")
    for subfolder_name in sorted(os.listdir(data_dir)):
        source = os.path.join(data_dir, subfolder_name)
        Compute_OC(source, subfolder_name, compute_OC, compute_top_ten, compute_repeatmatches, repeat_size)


