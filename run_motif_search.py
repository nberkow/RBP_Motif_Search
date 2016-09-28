import motifSearchRunner
import sys


if __name__ == "__main__":

	msr = motifSearchRunner.motifSearchRunner()

	#msr.set_login(sys.argv[1])
	#msr.set_pass(sys.argv[2])

	msr.set_background_from_fa('heptamers.fa', 'heptamers')
	msr.set_outpath(".")

	msr.run_motif_search(sys.argv[1:])

