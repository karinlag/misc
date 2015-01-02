import os.path
import sys
import glob
import re
from Bio import SeqIO

class MLST:
    def __init__(self, description):
        self.description = description
        self.mlst = self.description.split()[-1]
        self.mlstgene = self.mlst.split("_")[0]
        self.mlstno = self.mlst.split("_")[1]
        self.meta = "PERFECT"
        if "PERFECT" not in self.description:
            self.meta = "WARNING"

    def __str__(self):
        return "%s\t%s\t%s\t%s" % (self.mlst, self.mlstgene, self.mlstno, self.meta)


def processFsaFile(filename):
    MLSTs = {}
    for record in SeqIO.parse(filename, "fasta"):
        newMLST = MLST(record.description)
        MLSTs[newMLST.mlstgene] = newMLST
    return MLSTs


def processDirectory(inputDir):
    fsapattern = os.path.join(inputDir, "*genome.fsa")
    files = glob.glob(fsapattern)
    resultDict = {}
    for input in files:
        resultDict[input] = processFsaFile(input)
    return resultDict

def writeOutput(resultDict):
    keys = resultDict[resultDict.keys()[0]].keys()
    keys.sort()
    output = "SampleNo\tSampleName"
    for k in keys:
        output += "\t" + k + "\t" + k + "_NO" + "\t" + k + "_meta"
    print output
    for mlstFile in resultDict:
        mlsts = resultDict[mlstFile]
        name = os.path.basename(mlstFile)
        name = name.rstrip("_spades_scaffolds_Hit_in_genome.fsa")
        name = re.sub("Angen-bacDNA\d*-", "", name)
        name = name.split("_")[0]
        name = name.split("-")[0] + "\t" + "-".join(name.split("-")[1:])
        output = name
        for key in keys:
            output += "\t" + mlsts[key].mlst + "\t" + mlsts[key].mlstno + "\t" + mlsts[key].meta
        print output


if __name__=="__main__":
    inputDir = sys.argv[1]
    resultDict = processDirectory(inputDir)
    writeOutput(resultDict)