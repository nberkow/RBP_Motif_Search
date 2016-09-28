import motifSearcher
from Bio import SeqIO
import scipy
import scipy.stats as stat
import operator
import sys

class motifSearchRunner:


	def __init__(self):
		self.bg_set = False
		self.means_by_motif = {}
		self.sds_by_motif = {}
		self.n = 7 
		self.out_path = "."
		self.login = ''
		self.passwd = ''

	pnorm = stat.norm.cdf

	def set_outpath(self, op):
		self.out_path = op

	def set_login(self, login):
		self.login = login

	def set_pass(self, passwd):
		self.passwd = passwd

	def average_and_sd(self, distinct_nmers, motif_matches):


		# take into account the motif match scores and the number of times each
		# nmer occurred in the sequence set to find the mean and standard dev

		# Track all of these quantities for each motif. The motif_id is the key
		# for all dictionaries
		sums_by_motif = {}
		square_sums_by_motif = {}
		counts_by_motif = {}
		means_by_motif = {}
		sds_by_motif = {}
		z_scores_by_motif = {}

		for m in motif_matches.keys():

			match_scores = motif_matches[m]         # a list of scores, one for each motif with this nmer
			nmer_locations = distinct_nmers[m]      # all of the locations where the nmer appeared
			nmer_count = len(nmer_locations)        # the number of times the nmer appeared in the fasta			

			if match_scores:
				for match_score in match_scores:
					id = match_score.motifId
					if not sums_by_motif.has_key(id):
						sums_by_motif[id] = 0
						square_sums_by_motif[id] = 0
						counts_by_motif[id] = 0
						means_by_motif[id] = 0
						sds_by_motif[id] = 0
					counts_by_motif[id] = counts_by_motif[id] + nmer_count
					sums_by_motif[id] = sums_by_motif[id] + (match_score.getData()[1] * nmer_count)

		## Calculate the mean score for each motif
		for id in sums_by_motif.keys():
			sum = sums_by_motif[id]
			count = counts_by_motif[id]
			means_by_motif[id] = sum/count

		## Calculate the squares for the the standard deviations
		for m in motif_matches.keys():
			match_scores = motif_matches[m]
			nmer_locations = distinct_nmers[m]
			nmer_count = len(nmer_locations)
		
			if match_scores:
				for match_score in match_scores:
					id = match_score.motifId
					mean = means_by_motif[id]
					square_sums_by_motif[id] = square_sums_by_motif[id] + (((match_score.getData()[1] - mean)**2) * nmer_count)


		## Calculate the standard deviations
		for id in square_sums_by_motif.keys():
			square_sum = square_sums_by_motif[id]
			count = counts_by_motif[id]
			sds_by_motif[id] = (square_sum/(count - 1))**0.5

		## Calculate the z scores
		for m in motif_matches.keys():
			match_scores = motif_matches[m]
			if match_scores:
				for match_score in match_scores:
					id = match_score.motifId
					mean = means_by_motif[id]
					sd = sds_by_motif[id]
					score = match_score.getData()[1]
					z = (score - mean)/(sd**0.5)
					if not z_scores_by_motif.has_key(id):
						z_scores_by_motif[id] = []
					z_scores_by_motif[id].append(str(z))

		return (means_by_motif, sds_by_motif, z_scores_by_motif, counts_by_motif)

	def get_matches(self, distinct_nmers, species, n):
		mSearcher = motifSearcher.MotifSearch()
		sequence_matches = mSearcher.search_all_sequences(distinct_nmers, species, n)
		return sequence_matches

	def get_distinct_nmer_set(self, fasta_file, n):
		# get the set of distinct nmers. Record the sequence name and position
		# of each instance
		distinct_sequences = {}
		fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

		#print "getting distinct nmers"
		i = 0
		total = 0

    		for fasta in fasta_sequences:

			i += 1
			if i % 1000 == 0:
				pass
				#print i

        		name, sequence = fasta.id, fasta.seq.tostring()
			for x in range(len(sequence) - n + 1):
				total += 1
				nmer = sequence[x:(x+n)]
				if not distinct_sequences.has_key(nmer):
					distinct_sequences[nmer] = []
				distinct_sequences[nmer].append((name, total))

		print "%s total\n%s distinct" % (total, len(distinct_sequences.keys()))

		# return the set of distinct nmers and the total count of nmers checked
		return (distinct_sequences, total, i)

	def print_scores(self, distinct_nmers, motif_matches, out_file):
		print >> out_file, "motif ID\tmotif score\tposition\tnmer\tRBP Symbol\tSpecies\tnmer count"
		for m in distinct_nmers.keys():
			match_scores = motif_matches[m]
			nmer_count = len(distinct_nmers[m])
			if match_scores:
				for p in match_scores:
					print >> out_file, "%s\t%s" % (p, nmer_count)

	def call_significant_matches(self, distinct_nmers, motif_matches, means_by_motif, sds_by_motif, total_nmers, ap):
		# Test each distinct nmer for significance. Lookup the location information for the hits
		hits = []

		if ap == False:
			# Bonferonni correction for multiple testing
			alpha = 0.05
		
			total_matches = float(len(motif_matches.keys()) * (4**7))
			alpha_prime = alpha/total_matches

			print "total: %s" % (total_nmers)
			print "a1: %s" % (alpha_prime)
		else:
			alpha_prime = ap

		# Sidak correction for multiple testing
		#alpha = 0.05
		#total_matches = float(len(motif_matches.keys()) * total_nmers)
		#alpha_prime = 1 - (1 - alpha)**(1/total_matches)
		#print alpha_prime

		# call hits
		for nmer in motif_matches.keys():
			motif_match_set = motif_matches[nmer]
			if motif_match_set:
				for motif_match in motif_match_set:
					id = motif_match.motifId

					if id in means_by_motif.keys():	
						mean = means_by_motif[id]
						sd = sds_by_motif[id]
						n = len(motif_match_set)
						p_val = self.calculate_p_val(mean, sd, motif_match.motifScore, n)

						if p_val <= alpha_prime:
							locations = distinct_nmers[nmer]
							for loc in locations:
								(seq_name, pos) = loc
								(id, score, foo, seq, symb, species) = motif_match.getData() 
								hits.append((seq_name, pos, id, score, seq, symb, p_val))

		return hits

	def summarize_hits(self, hits, total_nmer_count):
		hit_counts = {}
		for h in hits:
			(seq_name, pos, id, score, seq, symb, p_val) = h
			if not hit_counts.has_key((id, symb)):
				hit_counts[(id, symb)] = 0
			hit_counts[(id, symb)] = hit_counts[(id, symb)] + 1
		sorted_hits = sorted(hit_counts.iteritems(), key=operator.itemgetter(1))
		sorted_hits.reverse()
		summary = "id\tsymbol\tmatches\tmatch frequency\n"
		for k in sorted_hits:
			id = k[0][0]
			symbol = k[0][1]
			count = k[1]
			norm_count = count/total_nmer_count
			summary = summary + "%s\t%s\t%s\t%s\n" % (id, symbol, count, float(count)/total_nmer_count)

		return summary 

	def get_seq_hit_counts(self, hits, total_seq_count):
		seq_report = "motif ID\tsymbol\tsequences matched\tfrequency\n"

		# keep track of which sequences have been counted for which motifs
		seq_sets_by_motif = {}

		# count the number of uniqe sequences that contain each motif
		seq_counts = {}

		for h in hits:
			(seq_name, pos, id, score, seq, symb, p_val) = h
			if not seq_sets_by_motif.has_key(id):
				seq_sets_by_motif[id] = set()

			if not seq_name in seq_sets_by_motif[id]:
				seq_sets_by_motif[id].add(seq_name)
				if not seq_counts.has_key((id, symb)):
					seq_counts[(id, symb)] = 0
				seq_counts[(id, symb)] += 1
							
		return seq_counts

	def get_seq_hit_count_summary(self, seq_counts, contrast_seq_counts, total_contrast_seq_count, total_seq_count):
		sorted_seqs = sorted(seq_counts.iteritems(), key=operator.itemgetter(1))
		sorted_seqs.reverse()

		seq_report = "motif ID\tsymbol\tsequences matched\tcontrast matched\tfrequency\tcontrast frequency\tenrichment\n"		

		for k in sorted_seqs:
			(id, symbol) = k[0]
			count = k[1]
			freq = (count + 0.)/total_seq_count

			if contrast_seq_counts.has_key((id, symbol)):
				contrast_count = contrast_seq_counts[(id, symbol)]
				contrast_freq = (contrast_count + 0.)/total_contrast_seq_count
				enrichment = freq/contrast_freq
			else:
				contrast_count = 0
				contrast_freq = 0.0
				enrichmrent = 'inf'
			
			seq_report = seq_report + "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (	id, symbol, count, 
											contrast_count, freq, contrast_freq, 
											enrichment)
		return seq_report



	def format_hits(self, hits):
		hit_report = "seq name\tposition\tID\tscore\tsequence\tsymbol\tp_val\n"
		for h in hits:
			hit_report += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % h 
		return hit_report
	
	def calculate_p_val(self, mean, sd, score, n):
		z = (mean-score)/(sd**0.5)
		p_val = self.pnorm(-abs(z))
		return p_val

	def set_background_from_fa(self, fa_file, bg_prefix):
		
		## get the set of distinct nmmer and the total number of nmers in the background fasta
		(distinct_nmers, total_nmer_count, total_bg_seq_count) = self.get_distinct_nmer_set(fa_file, self.n)

	
		## get all of the pairwise matches between the background nmer set and the motif set
		species = "('Homo_sapiens', 'Danio_rerio', 'Mus_musculus', 'Arabidopsis_thaliana')"
		#species = "('Arabidopsis_thaliana', 'Physcomitrella_patens')"

		motif_matches = self.get_matches(distinct_nmers, species, self.n)
		## get the average match score and stdev for each motif based on the background
		(means_by_motif, sds_by_motif, z_scores_by_motif, counts_by_motif) = self.average_and_sd(distinct_nmers, motif_matches)

		outfile = open("%s_background_scores_bnf.txt" % (bg_prefix), 'w')
		self.print_scores(distinct_nmers, motif_matches, outfile)
		outfile.close()

		## report the background statistics
		outfile = open("%s_background_stats_bnf.txt" % (bg_prefix), 'w')
		print >> outfile, "motif\tmean\tsd"
		for m in means_by_motif:
			mean = means_by_motif[m]
			sd = sds_by_motif[m]
			print >> outfile, m, mean, sd

		outfile = open("%s_background_zscores_bnf.txt" % (bg_prefix), 'w')
		for k in z_scores_by_motif.keys():
			z = z_scores_by_motif[k]
			print >> outfile, k + "," + ",".join(z)

		self.means_by_motif = means_by_motif
		self.sds_by_motif = sds_by_motif
		self.bg_set = True
		outfile.close()

	def set_background_from_summary(self, summary_file_name, bg_prefix):
		sf = open(summary_file_name, 'r')
		sf.readline()
		for s in sf:
			(id, mean, sd) = s.split(" ")	
			self.means_by_motif[id] = float(mean)
			self.sds_by_motif[id] = float(sd)	

		sf.close()
		self.bg_set = True


	def run_motif_search(self, args):
		## Find the background distribution from a fasta file
		self.n = 7
		target_files = args

		if self.bg_set:

			## Call significant scores in each test set. Use the first file as a contrast for the rest

			contrast_hits = None
			contrast_seq_count = 0
			for f_name in target_files:
				target_prefix = "".join(f_name.split("/")[-1].split(".")[0:-1])

				print "getting distinct nmer set"
				(distinct_nmers, total_nmer_count, total_seq_count) = self.get_distinct_nmer_set(f_name, self.n)

				#species = "('Homo_sapiens')"
				print "finding motif matches"
				#species = "('Homo_sapiens', 'Danio_rerio', 'Mus_musculus')"
				species = "('Arabidopsis_thaliana', 'Physcomitrella_patens')"

				motif_matches = self.get_matches(distinct_nmers, species, self.n)

				print "calling hits"
				hits = self.call_significant_matches(distinct_nmers, motif_matches, self.means_by_motif, 
								self.sds_by_motif, total_nmer_count, False)

				print "processiong contrast set"
				if not contrast_hits:
					contrast_hits = hits
					total_contrast_seq_count = total_seq_count
					total_contrast_nmer_count = total_nmer_count
					contrast_seq_counts = self.get_seq_hit_counts(contrast_hits, total_contrast_seq_count)

				print "reporting hits"
				hits_file = open("%s/%s.hits_bnf.txt" % (self.out_path, target_prefix), "w")
				hits_report = self.format_hits(hits)
				print >> hits_file, hits_report

				print "summarizing hits"
				summary = self.summarize_hits(hits, total_nmer_count)
				summary_file = open("%s/%s.summary_bnf.txt" % (self.out_path, target_prefix), "w")
				print >> summary_file, summary

				print "preparing sequence summary"
				seq_counts = self.get_seq_hit_counts(hits, total_seq_count)
				seq_summary = self.get_seq_hit_count_summary(seq_counts, contrast_seq_counts, total_seq_count, total_contrast_seq_count)
				seq_summary_file = open("%s/%s.seq_summary_bnf.txt" % (self.out_path, target_prefix), "w")
				print >> seq_summary_file, seq_summary

				summary_file.close()
				hits_file.close()
				seq_summary_file.close()

		else:
			print "no background loaded. exiting"















