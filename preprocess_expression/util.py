import sys
import string
import subprocess
import os


DNA_COMP = None

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global DNA_COMP

    if DNA_COMP is None:
        DNA_COMP = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(DNA_COMP)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]

        
def sort_bam(input_bam, output_prefix):
    """Calls samtools sort on input_bam filename and writes to
    output_bam. Takes into account that the command line arguments 
    for samtools sort have changed between versions."""

    output_bam = output_prefix + ".sort.bam"
    
    # first try new way of using samtools sort
    failed = False
    cmd = "samtools sort -o " + output_bam + " " + input_bam
    sys.stderr.write("running command: %s\n" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        sys.stderr.write("samtools sort command failed:\n%s\n" %
                         str(e))
        failed = True
    if not os.path.exists(output_bam):
        sys.stderr.write("output file %s does not exist\n" % output_bam)
        failed = True
        
    if failed:
        # OLD way of calling samtools (changed in newer versions)
        sys.stderr.write("samtools sort command failed, trying old samtools "
                         "syntax\n")
        
        cmd = "samtools sort " + input_bam + " " + output_prefix
        sys.stderr.write("running command: %s\n" % cmd)

        try:
            subprocess.check_call(cmd, shell=True)
        except Exception as e:
            sys.stderr.write("samtools sort command failed:\n%s\n" %
                             str(e))
            exit(1)
        
        if not os.path.exists(paths.sorted_output_bam):
            raise IOError("Failed to create sorted BAM file '%s'" %
                          paths.sorted_output_bam)



def is_gzipped(filename):
    """Checks first two bytes of provided filename and looks for
    gzip magic number. Returns true if it is a gzipped file"""
    f = open(filename, "rb")

    # read first two bytes
    byte1 = f.read(1)
    byte2 = f.read(1)
    
    f.close()

    # check against gzip magic number 1f8b
    return (byte1 == chr(0x1f)) and (byte2 == chr(0x8b))



def check_pysam_version(min_pysam_ver="0.8.4"):
    """Checks that the imported version of pysam is greater than
    or equal to provided version. Returns 0 if version is high enough,
    raises ImportWarning otherwise."""
    import pysam

    min_ver = [int(x) for x in min_pysam_ver.split(".")]
    pysam_ver = [int(x) for x in pysam.__version__.split(".")]

    n_ver = min(len(pysam_ver), len(min_pysam_ver))
    
    for i in range(n_ver):
        if pysam_ver[i] < min_ver[i]:
            raise ImportWarning("pysam version is %s, but pysam version %s "
                                "or greater is required" % (pysam.__version__,
                                min_pysam_ver))
        if pysam_ver[i] > min_ver[i]:
            # version like 1.0 beats version like 0.8
            break
        
    return 0
        
    
                                

    
