import sys

#########################################################################################################################################################################
# Script created by Jens Theine, 2020
# Script to renumber all contigs (identified by ">"-character stating a new line).
# To use the script use f.ex.: "python Renametigs.py --in /home/juser/assembly.fasta --out /home/juser/assembly_reordered_tigs.fasta --prefix fabulous_contig"
#########################################################################################################################################################################

__usage__ = """
        python renametigs.py
        --in <FULL/PATH/TO/ASSEMBLY/FILE.fasta>
        --out <FULL/PATH/TO/OUTPUT/FILE.fasta>
        --prefix <CONTIG-PREFIX>
        """

def main(parameters):
        input = parameters[parameters.index('--in')+1]
        output = parameters[parameters.index('--out')+1]
        prefix = parameters[parameters.index('--prefix')+1]

        i = 0

        open(output, 'w').close()

        with open(input, 'r') as f:
                for line in f:
                        if ">" in line:
                                #oldcontig = line.split('>')[1]
                                i = i + 1
                                contig = ">" + prefix + '{0:04}'.format(i)
                                #print contig
                                with open(output, "a") as lineout:
                                        lineout.write(contig + "\n")
                        else:
                                #print line
                                with open(output, "a") as lineout:
                                        lineout.write(line)

        print "\ndone."
        print "- changed names in " + str(i) + " contigs and saved to " + output + "."

if __name__ == '__main__':
        if '--in' in sys.argv and '--out' in sys.argv and '--prefix' in sys.argv:
                main( sys.argv )
        else:   
                sys.exit( __usage__ )
