#!/usr/bin/env python
#-------------------------------------------------------------
# Author: Thomas Schwarzl <schwarzl@embl.de>
# With the help of: Christian Hauer <chauer@embl.de>
# Licenced under MIT Creative Licence
# Last change: 21 October 2014
#-------------------------------------------------------------

#-------------------------------------------------------------
# pip install numpy --user
import sys, os, re, csv, random, math, logging, shutil, gzip, traceback
import datetime, time

try:
    import HTSeq
except Exception:
    print "Please install the HTSeq framework e.g. like this"
    print "pip install HTSeq"
    print "pip install HTSeq --user"
    os._exit(1)

from pybedtools import BedTool

from sortedcontainers import SortedSet
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from random import randint   # can be removed in the real script
from itertools import islice



#-------------------------------------------------------------
class iCLIP:
   data          = {}
   features      = {}
   reads         = {} #probably to delete
   stats         = {}
   sites_unique  = {}
   sites         = {}
   stranded      = True
   bamfile       = ""
   gtfname       = ""
   data          = {}
   datadf        = {}
   verbose       = False
   deletionRatio = 0
   minDeletions  = 1
   halfwindow    = 100
   sortOutput    = True
   primary       = False
   filterExons   = True
   maxReadIntervalLength = 0
   minAlignmentQuality = 0
   
   maxReadLength = 0
   minReadLength = 0
   sites_unique  = 0

   files         = {}
   readLength    = {}

   coverage      = {}

   # Constructor
   # Reads all the options and creates the data structures needed 
   def __init__(self, options):
      if hasattr(options, 'bamfile'):
           self.bamfile                  = options.bamfile

      self.outputDir                = options.outputDir

      self.bedOutputDir             = os.path.join(options.outputDir, 'bed')

      if hasattr(options, 'maxReadIntervalLength'):
        self.maxReadIntervalLength    = options.maxReadIntervalLength
      if hasattr(options, 'minAlignmentQuality'):
        self.minAlignmentQuality      = options.minAlignmentQuality

      self.stranded                 = options.stranded

      self.verbose                  = options.verbose

      if hasattr(options, 'deletionRatio'):
        self.deletionRatio            = options.deletionRatio

      if hasattr(options, 'minDeletions'):
        self.minDeletions             = options.minDeletions

      if hasattr(options, 'gtfname'):
        self.gtfname                  = options.gtfname

      if hasattr(options, 'halfwindow'):
        self.halfwindow               = options.halfwindow

      if hasattr(options, 'sortOutput'):
        self.sortOutput               = options.sortOutput

      if hasattr(options, 'primary'):
        self.primary                  = options.primary

      if hasattr(options, 'filterExons'):
         self.filterExons              = options.filterExons

      self.features  =  { 'exons':  HTSeq.GenomicArrayOfSets( "auto", stranded=self.stranded ) ,
                        'introns':  HTSeq.GenomicArrayOfSets( "auto", stranded=self.stranded ) }


      #self.reads = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)

      self.stats =  {'total_reads':                 0,
                     'above quality criteria':      0,
                     'reads with deletions':        0,
                     'reads with mismatches':       0,
                     'reads with insertions':       0,
                     'total deletions':             0,
                     'total insertions':            0,
                     'total mismatches':            0,
                     'average deletion length':     0,
                     'average insertion length':    0,
                     'average mismatch length':     0,
                     'maxReadInterval':             self.maxReadIntervalLength,
                     'minAlignmentQuality':         self.minAlignmentQuality}


      # create the output dirs
      self.createDir(self.outputDir)
      self.createDir(self.bedOutputDir)



   # =====================================================================================
   # Destructor closes the file connections
   def __exit__(self, type, value, traceback):
      for key, file in self.files:
         file.close()
      
   # =====================================================================================
   def run(self):   
      # calculateDistancesFromReads
      self.calculateDistancesFromReads()
      self.calculateDistancesFromSites()

   # =====================================================================================
   # The Read lengths will be determined, the starting positions
   # of all reads and the positions of deletions
   def readBamFile(self):

      # initiate file handlers for the deletion and insertion bed
      self.files['deletion']  = gzip.open(os.path.join(self.bedOutputDir, 'cims-deletions.bed.gz'), 'w')
      self.files['insertion'] = gzip.open(os.path.join(self.bedOutputDir, 'cims-insertions.bed.gz'), 'w')
      self.files['reads']     = gzip.open(os.path.join(self.bedOutputDir, 'cims-reads.bed.gz'), 'w')


      # get a BAM Reader from HTSeq
      almnt_file = HTSeq.BAM_Reader( self.bamfile )
      
      print "Reading the Bam File"

      # Read through the BAM file and retrieve aligned reads
      for almnt in almnt_file:
            # Counts the read
            self.stats['total_reads'] += 1

            if self.stats['total_reads'] % 10000 == 0:
                print "Read in " + str(self.stats['total_reads']) + " reads"

            if self.readFullfillsQualityCriteria(almnt):
                      # Counts the read if it is above a certain quality criteria
                      self.stats['above quality criteria'] += 1

                      # Calculate read length and also the minimum and maximum read length
                      readLength = self.calcMinMax(almnt)

                      # Calculate the adapter length and add it to a table: adapter length versus read length
                      self.increaseDataCount('adapters', self.getAdapterSequence(almnt.read.name), readLength)

                      # Calculate the starting position of the read and store it in self.sites
                      self.increaseSiteCount('start', self.determineStartSite(almnt.iv), almnt.read.name)

                      # Calculate the middle position of the read and store it in self.sites
                      self.increaseSiteCount('middle', self.determineMiddleSite(almnt.iv), almnt.read.name)

                      # Calculate the end position of the read and store it in self.sites
                      self.increaseSiteCount('end', self.determineEndSite(almnt.iv), almnt.read.name)

                      # Calculate the position of insertions, deletion and mismatches and store them in self.sites
                      self.addCigarInformation(almnt, readLength)

                      # Calculate the cross link site - this will be in future either the start site
                      # or for example a deletion site if available
                      # self.increaseSiteCount('crosslink', self.determineCrosslinkSite(almnt.iv), almnt.read.name)

                      # Create a coverage plot and store it in self.sites
                      # self.increaseSiteCount('coverage', almnt.iv, almnt.read.name)

                      # Store the individual reads in self.reads
                      #self.addRead(almnt, readLength)

                      # Stores the sequence Frequency per base
                      self.countSequenceFrequencies(almnt.read.seq, readLength)

                      # Write the Read as CIMS bed file
                      self.writeReadToCIMSBed(almnt)

                      # store the read length
                      self.readLength[almnt.read.name] = readLength

                      # stores the read counts for read length
                      self.increaseDataCount('readcount', readLength, "counts")

      # Do some stats calculations after all reads have been read.
      self.calculateReadStats()


   #-------------------------------------------------------------
   # This function calculates the distances of a read
   # to certain features. Therefore it will run through all the
   # reads of the BAM file again. This could be also stored in
   # the memory if that would be faster and memory not an issue
   def calculateDistancesFromReads(self):
      # Looping through the reads again and calculating the distances
      print "Calculating Distances From Reads"

      # get a HTSeq bam file reader
      almnt_file = HTSeq.BAM_Reader( self.bamfile )

      i = 1
      # run through the aligned reads
      for almnt in almnt_file:

         i += 1
         if i % 10000 == 0:
             print "%s reads processed" % i

         # if the criteria is fullfilled
         if self.readFullfillsQualityCriteria(almnt):

            # calculate the features of the read
            readLength  = self.getSequenceLength(almnt)
            startSite    = self.determineStartSite(almnt.iv)
            middleSite   = self.determineMiddleSite(almnt.iv)
            endSite      = self.determineEndSite(almnt.iv)

            # calculate the distance of the read to start certain features.
            self.doCalculateDistance(startSite,  'start',      readLength, 'start-to-start')
            self.doCalculateDistance(startSite,  'deletions',  readLength, 'start-to-deletion')
            self.doCalculateDistance(middleSite, 'deletions',  readLength, 'middle-to-deletion')
            self.doCalculateDistance(endSite,    'deletions',  readLength, 'end-to-deletion')
            self.doCalculateDistance(startSite,  'insertions', readLength, 'start-to-insertion')
            self.doCalculateDistance(middleSite, 'insertions', readLength, 'middle-to-insertion')
            self.doCalculateDistance(endSite,    'insertions', readLength, 'end-to-insertion')

   #-------------------------------------------------------------
   # This function calculates the distance from features to other
   # features. So e.g. from one deletion site to another deletion
   # site.
   def calculateDistancesFromSites(self):
      print "Calculating Distances From Sites"


      #for (start, value) in self.sites['start'].steps():
      #   if len(value) > 0:
      #      self.doCalculateDistance(start,  'start',          value, 'start-to-start')

      for (exonIntron, value) in self.sites['exon-intron'].steps():
         if len(value) > 0:
            self.doCalculateDistance(exonIntron,  'deletions',      value, 'exon-intron-to-deletion')
            self.doCalculateDistance(exonIntron,  'insertions',     value, 'exon-intron-to-insertion')
            self.doCalculateDistance(exonIntron,  'start',          value, 'exon-intron-to-start')
            self.doCalculateDistance(exonIntron,  'end',            value, 'exon-intron-to-end')
            self.doCalculateDistance(exonIntron,  'middle',         value, 'exon-intron-to-middle')
            self.doCalculateCoverage(exonIntron,  value, 'exon-intron-coverage')

      for (intronExon, value) in self.sites['intron-exon'].steps():
         if len(value) > 0:
            self.doCalculateDistance(intronExon,  'deletions',      value, 'intron-exon-to-deletion')
            self.doCalculateDistance(intronExon,  'insertions',     value, 'intron-exon-to-insertion')
            self.doCalculateDistance(intronExon,  'start',          value, 'intron-exon-to-start')
            self.doCalculateDistance(intronExon,  'end',            value, 'intron-exon-to-end')
            self.doCalculateDistance(intronExon,  'middle',         value, 'intron-exon-to-middle')
            self.doCalculateCoverage(intronExon,  value, 'intron-exon-coverage')

      for (geneStart, value) in self.sites['gene-start'].steps():
         if len(value) > 0:
            self.doCalculateDistance(geneStart,  'deletions',       value, 'gene-start-to-deletion')
            self.doCalculateDistance(geneStart,  'insertions',      value, 'gene-start-to-insertion')


      for (geneEnd, value) in self.sites['gene-end'].steps():
         if len(value) > 0:
            self.doCalculateDistance(geneEnd,  'deletions',         value, 'gene-end-to-deletion')
            self.doCalculateDistance(geneEnd,  'insertions',        value, 'gene-end-to-insertion')

      for (deletion, value) in self.sites['deletions'].steps():
         if len(value) > 0:
             self.doCalculateDistance(deletion,  'deletions',       value, 'deletion-to-deletions')
             self.doCalculateDistance(deletion,  'insertions',      value, 'deletion-to-insertions')
             self.doCalculateDistance(deletion,  'start',           value, 'deletion-to-start')
             self.doCalculateDistance(deletion,  'middle',          value, 'deletion-to-middle')
             self.doCalculateDistance(deletion,  'end',             value, 'deletion-to-end')

      for (insertion, value) in self.sites['insertions'].steps():
         if len(value) > 0:
            self.doCalculateDistance(insertion,  'insertions',      value, 'insertion-to-insertions')

      # calculate the ratios
      #self.doCalculateDeletionRatio('exon-intron-to-deletion', 'exon-intron-coverage', offset=1)
      #self.doCalculateDeletionRatio('exon-intron-to-deletion', 'exon-intron-coverage', offset=1)


   # #-------------------------------------------------------------
   # # calculate the deletion rate per base
   # def doCalculateDeletionRatio(self, mutationKey, coverageKey, offset=0):
   #      for dist, readLengthDict in self.data[mutationKey]:
   #          for readLength, counts in readLengthDict:
   #              val = offset
   #              try:
   #                  val = self.data[coverageKey][dist][readLength]
   #              except IndexError:
   #                  val = offset
   #
   #              if val <= 0:
   #                  logging.warning("Offset ")
   #              int(dist)
   #
   #             # TODO 1) finish the offset calculation
   #                    2) check if the distances are calculated correctly
   #
   #
   #          # deletions are not counted in the coverage!
   #
   #
   #      for i in self.data['exon-intron-to-deletion']:
   #          pass

   #-------------------------------------------------------------
   # This function calculates Distances from a site
   # site: HTSeq.GenomicPosition of the site
   # sitekey: look for sites for distance calculation in that specific container
   # value: ?
   # datakey: name of the data storage for this distance

   def doCalculateDistance(self, site, sitekey,  value, datakey, min=False):
      # print "doCalculateDistance %s %s %s %s " % (str(site), sitekey, marker, datakey)

      # if min == True it will only return one item
      # returns list of tuples: distance to readLength
      distList = self.calcDistance(sitekey, site, value, min)

      if len(distList) == 0:
          self.increaseDataCount(datakey, "unknown", "unknown")
      else:
          for dist, readLength in distList:
               #print "Storing %s and with read length %s" % (dist, readLength)
               self.increaseDataCount(datakey, str(dist), readLength)


   #-------------------------------------------------------------
   # get the background coverage for the window
   def doCalculateCoverage(self, eventSite,  eventValue, datakey):
      window = HTSeq.GenomicInterval(eventSite.chrom,
                                     max(0, eventSite.start - self.halfwindow),
                                     eventSite.end + self.halfwindow,
                                     eventSite.strand)

# TODO it does not work
      #   File "/g/hentze/projects/iCLIP/BindingSites/src/iclipper/lib/iCLIP.py", line 346, in doCalculateCoverage
 #   for interval, value in cov[window].steps():
# TypeError: string indices must be integers, not HTSeq._HTSeq.GenomicInterval

      # # returns genomic interval (interval) and a set (values)
      # for readLength, cov in self.coverage:
      #     for interval, value in cov[window].steps():
      #         # get the distance from the start to the middle of the event
      #         dist = math.floor(interval.start - (eventSite.start + eventSite.end) / 2)
      #
      #         # add the interval to the data table
      #         for i in range(0, interval.length - 1):
      #             self.increaseDataCount(datakey, dist + i, readLength, value)

   #-------------------------------------------------------------
   # Stores the nucleotide frequencies per read length
   def countSequenceFrequencies(self, sequence, readLength):
      pos = 1
      for nucleotide in sequence:
         self.increaseDataCount("nucleotides-" + str(readLength), nucleotide, str(pos))
         self.increaseDataCount("nucleotides-all", nucleotide, str(pos))
         pos += 1

   #-------------------------------------------------------------
   def calcDistance(self, sitekey, eventSite, eventValue, min=False):
      retList = []

      #print "[[[ %s, %s, %s ]]]" % (sitekey, eventSite, eventValue)
      # Get the window for the specific site and get all features in there as GenomicArray
      # Then we iterate over this Genomic array which gives back an interval and a value
      # The value corresponds to how many of these events occur at the site
      #print "Strand: %s  eventsite %s eventvalue %s" % (eventSite., eventSite.end, eventSite, eventValue)
      window = HTSeq.GenomicInterval(eventSite.chrom,
                                                max(0, eventSite.start_d - self.halfwindow),
                                                eventSite.end_d + self.halfwindow,
                                                eventSite.strand)

      defaultValue =  self.halfwindow + 1

      # iterable
      i = defaultValue

      # returns genomic interval (interval) and a set (values)
      for interval, values in self.sites[sitekey][window].steps():
         # the sets can be empty, therefor check and
         # go through every taken site
         if len(values) > 1:
            # strings are returned when there is only one value, avoid this
            # by creating an iterable element (list) manually
            if type(values) is str:
                values = [ values ]

            #print "VALUES: %s " % values
            # now iterate through all the sites
            for val in values:
                # calculate the distance for each value which is not on the same read.
                if val != eventValue:
                    # get the read length
                    readLength = self.readLength[val]

                    # calculate distance
                    dist = math.floor((interval.start_d + interval.end_d) / 2 - (eventSite.start_d + eventSite.end_d) / 2)
                    #print " === distance %s" % dist

                    # if the minimum is seeked, just look for the first unique one, since everything is ordered anyway
                    if min:
                        # if the stored is smaller than the current one, then jump out of the loop and return the stored
                        # value
                        if abs(i) > dist:
                            return [(dist, readLength)]
                        # else store this values as new value
                        else:
                            i = dist
                    else:
                        retList.append((dist, readLength))


      # if there is no read greater downstream, then take the best upstream
      if min and abs(i) < defaultValue:
          retList.append((i, readLength))

      return retList
   
   #-------------------------------------------------------------
   # Gets the steps for an interval in a site GenomicArray
   # which can be indexed by sitekey
   # the start is set so it cannot be a negativ number
   def getSitesInWindow(self, sitekey, startSite, halfwindow):
      #print "%s Start %s, Site: %s, Stop %s" % (startSite.strand, max(0, startSite.start - halfwindow), startSite.start, startSite.end + halfwindow)
      return self.sites[sitekey][ HTSeq.GenomicInterval(startSite.chrom, max(0, startSite.start - halfwindow), startSite.end + halfwindow, startSite.strand) ]
      
   #-------------------------------------------------------------
   # Filter Deletions:
   # Adds a GenomicArray to sites under key 'selected-deletions'
   # with a minimum of 'options.minDeletions' reads and
   # a ratio of 'options.deletionRatio' deletions to
   # background
   def filterDeletions(self):
         deletionRatio = self.deletionRatio
         
         for iv, value in self.sites['deletions'].steps():
            if value >= self.minDeletions:
               #print "----"
               #print sites['coverage'][iv]
               # TODO: Select for ratio
               self.sites['selected-deletions'][iv] += value

      
   #-------------------------------------------------------------
   # Parse Cigar string and add to the statistics
   def addCigarInformation(self, almnt, readLength):
       variations = self.parseCigar(almnt)
       
       if len(variations['deletions']) > 0:
            self.increaseStats('reads with deletions')
       if len(variations['deletions']) > 0:
            self.increaseStats('reads with insertions')
             
       for variation in variations['deletions']:
            self.increaseStats('total deletions')
            self.increaseStats('average deletion length', variation.ref_iv.length)
            self.increaseSiteCount('deletions', variation.ref_iv, almnt.read.name)
            self.increaseDataCount('deletions', variation.ref_iv.length, readLength)
            self.addPositionsOnRead('deletion', variation.query_from, variation.query_to + 1, readLength, almnt)

       for variation in variations['insertions']:
            self.increaseStats('total insertions')
            self.increaseStats('average insertion length', variation.ref_iv.length)
            self.increaseSiteCount('insertions', self.getGenomicIntervalWithEndOffset(variation.ref_iv, 1), almnt.read.name)
            self.increaseDataCount('insertions', variation.ref_iv.length, readLength)
            self.addPositionsOnRead('insertion', variation.query_from, variation.query_to, readLength, almnt)

       for variation in variations['hits']:
            # add to coverage plot
            self.addCoverage(variation.ref_iv, readLength)

   #-------------------------------------------------------------
   def getGenomicIntervalWithEndOffset(self, interval, offset):
       if interval.strand == "-":
           return(HTSeq.GenomicInterval(interval.chrom, interval.start, interval.end + 1, interval.strand))
       else:
           return(HTSeq.GenomicInterval(interval.chrom, interval.start - 1, interval.end, interval.strand))

   #-------------------------------------------------------------
   def addPositionsOnRead(self, key, query_from, query_to, readLength, almnt):
        self.writeVariationToCIMSBed(key, query_from, query_to, almnt)

        for i in range(query_from, query_to):
            self.increaseDataCount(str(key) + "-sites-perbase", almnt.read.seq[i-1], readLength)
            self.increaseDataCount(str(key) + "-sites",         i,                   readLength)

   #-------------------------------------------------------------
   # This function write a Bed file need for the CIMS/CITS analysis
   def writeVariationToCIMSBed(self, key, query_from, query_to, almnt):
       self.files[key].write(str(almnt.iv.chrom) + "\t" + str(almnt.iv.start) + "\t" + str(almnt.iv.end)  + "\t" + str(almnt.read.name) + "\t" + str(query_from) + "\t" + str(almnt.iv.strand) + "\n")

   #-------------------------------------------------------------
   # This function write a Bed file need for the CIMS/CITS analysis
   def writeReadToCIMSBed(self, almnt):
       self.files['reads'].write(str(almnt.iv.chrom) + "\t" + str(almnt.iv.start) + "\t" + str(almnt.iv.end) + "\t" + str(almnt.read.name) + "\t" + str(1) + "\t" + str(almnt.iv.strand) + "\n")

   #-------------------------------------------------------------
   # This method determines if a read fullfills the critera to be included in the analysis
   def readFullfillsQualityCriteria(self, almnt):
      if almnt.paired_end and almnt.pe_which == "second":
          return False
      else:
          return(almnt.aligned and
                 almnt.iv.length <= self.maxReadIntervalLength and
                 almnt.aQual >= self.minAlignmentQuality and
                 not almnt.failed_platform_qc and # SAM flag 0x0200
                 self.primaryFilter(almnt))

   #-------------------------------------------------------------
   # If the primary filter is activated (option primary is set),
   # for multimapping reads it will filter out the best location
   # using the SAM flag 0x0100
   def primaryFilter(self, almnt):
        if self.primary:
            return(almnt.not_primary_alignment)
        else:
            return(True)

   # ------------------------------------------------------------
   # Print verbose.
   # Prints the message
   def printv(self, string):
      if self.verbose == True:
         print "%s" % string
   

   #-------------------------------------------------------------
   def increaseStats(self, key, value=1):
        self.stats[key] += value
   
   #-------------------------------------------------------------
   # increases the data count by value
   def increaseDataCount(self, key, x, y, value=1):
        if key not in self.data:
            self.data[key] = { x: { y: value }}
        else:
            try:
                self.data[key][x][y] += value
            except KeyError:
                try:
                    self.data[key][x][y] = value
                except KeyError:
                    self.data[key][x] = { y: value }

   #-------------------------------------------------------------
   # increases the count of a site type by a given value
   def increaseSiteCount(self, key, pos, value=1):
      s = self.sites.get(key, HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)) #, storage="memmap")

      if (pos.end - pos.start) <= 0:
          logging.warning("Wanted to add site count " + key + " for " + pos.chrom + ":" + str(pos.end) + " - " + str(pos.start) + " with value " + str(value))
      else:
          try:
              s[pos] = value
              self.sites[key] = s
          except IndexError:
              logging.warning("IndexError when adding " + str(pos) + " value " + str(value) + " to key " + str(key))

   #-------------------------------------------------------------
   # increases the count of a site type by a given value
   def addCoverage(self, interval, readLength, value=1):
      readLength = str(readLength)
      s = self.coverage.get(readLength, HTSeq.GenomicArray("auto", stranded=self.stranded, typecode='i')) #, storage="memmap")

      try:
          #print "Adding %s  with value %s " % (interval, value)
          s[interval] += value
          self.coverage[readLength] = s
      except IndexError, KeyError:
          logging.warning("IndexError when adding coverage " + str(interval) + " value " + str(value) + " to key " + readLength)


   # #-------------------------------------------------------------
   # # adds a read
   # def addCoverage(self, almnt):
   #    if (almnt.iv.end - almnt.iv.start) > 0:
   #       self.reads[almnt.iv] += almnt.read.name
   
   #-------------------------------------------------------------
   # returns a list of GenomeIntervals for the positions of
   # deletions.
   def parseCigar(self, almnt):
      variations = { 'deletions': list(),
                     'insertions': list(),
                     'hits': list() }
      
      for i in almnt.cigar:
            if i.type == 'D':
                  variations['deletions'].append(i)
            elif i.type == 'I':
                  variations['insertions'].append(i)
            elif i.type == 'M':
                  variations['hits'].append(i)
         
      return(variations)

   #-------------------------------------------------------------
   # returns GenomicPosition for start site
   def determineStartSite(self, iv):
         return(HTSeq.GenomicPosition(iv.chrom, iv.start_d, iv.strand))
         
   #-------------------------------------------------------------
   # returns GenomicPosition for end site
   def determineEndSite(self, iv):
         return(HTSeq.GenomicPosition(iv.chrom, iv.end_d, iv.strand))
   
   #-------------------------------------------------------------
   # returns GenomicPosition for middle site
   # TODO: Check how this performs with gapped aligments
   def determineMiddleSite(self, iv):
         return(HTSeq.GenomicPosition(iv.chrom, round((iv.start_d + iv.end_d) / 2), iv.strand))

   #-------------------------------------------------------------
   # TODO: determine the crosslink site, means a basepair down of the start site
   def determineCrosslinkSite(self, iv):
      pos = 0
      
      if(self.stranded):
            if(iv.strand == "+"):
               pos = iv.start_d - 1
            elif(iv.strand == "-"):
               pos = iv.start + 1
            else:
               raise("Strand not known %s" % iv.strand)
      else:
            pos = iv.start_d - 2
      
      return(HTSeq.GenomicPosition(iv.chrom, pos, iv.strand))


   #-------------------------------------------------------------
   # this will convert the data structure to an output data
   # structure
   def convertDataToDataFrame(self):
       for i in self.data:
           self.datadf[i] = pd.DataFrame(self.data[i]).T.fillna(0)

   #-------------------------------------------------------------
   # sorts the data frame so the output will look nice
   def sortData(self):
      for key in self.datadf:
            self.datadf[key] = self.datadf[key].sort_index(axis=0, ascending=True)
            self.datadf[key] = self.datadf[key].sort_index(axis=1, ascending=True)

   #-------------------------------------------------------------
   # Write the Results to OutputDirectory and copies the
   # R script there, as well as executes it
   def writeOutput(self):
      print "Writing Output:"
      print "# Converting data to matrices and sorting"
      self.convertDataToDataFrame()
      if self.sortOutput:
         self.sortData()

      print "# Writing data tables"
      for keys in self.datadf:
         self.writeDataToFile(keys)
         
      print "# Writing the stats"
      self.writeStatsToFile(self.outputDir + '/stats.csv')

      print "# Writing Bed Graphs"
      for keys in self.sites:
         self.writeBedGraph(keys)

      print "# Copying the R file into the folder"
      # TODO
      # rfile = os.path.join('lib', 'iclipper.R')
      rfile = '/g/hentze/projects/iCLIP/BindingSites/src/iclipper/iclipper.R'

      if os.path.exists(rfile):
      #assert not os.path.isabs(rfile)
        dstdir = os.path.join(self.outputDir, 'iclipper.R')
        shutil.copy(rfile, dstdir)


   #-------------------------------------------------------------
   def writeStatsToFile(self, statsname):
         w = csv.writer(open(statsname, "w"))
         for key, val in self.stats.items():
               w.writerow([key, str(val)])

   #-------------------------------------------------------------
   def getSequenceLength(self, almnt):
      return(len(almnt.read.seq))
      
   #-------------------------------------------------------------
   # calculates the read length and stores the min and max length
   def calcMinMax(self, almnt):
      length = self.getSequenceLength(almnt)
      
      if not (self.maxReadLength == 0 and self.minReadLength == 0):
            if length > self.maxReadLength:
               self.maxReadLength = length
            if length < self.minReadLength:
               self.minReadLength = length
      return(length)

   #-------------------------------------------------------------
   # Extract the adapter sequence of a read name.
   # Here, we presume that it is the character string at the end of the name
   def getAdapterSequence(self, readname ):
      adapter = ""
      m      = re.search('[A-Z]+$', readname)
      if m:
         return(m.group(0))
      else:
         return("unknown")

   #-------------------------------------------------------------
   # Calculate read statistics
   def calculateReadStats(self):
      if self.stats['total deletions'] > 0:
         self.stats['average deletion length']  = self.stats['average deletion length']   / self.stats['total deletions']
      if self.stats['total insertions'] > 0:
         self.stats['average insertion length'] = self.stats['average insertion length']  / self.stats['total insertions']


   #-------------------------------------------------------------
   # Write the table from data to the output directory
   def writeDataToFile(self, slot):
         file = open(self.outputDir + "/data_" + slot + ".txt", "w")
         file.write(self.datadf[slot].to_string())
         file.close()

   #-------------------------------------------------------------
   # Write Bed Graph
   def writeBedGraph(self, index):
      if self.stranded:
            self.getGenomicArrayFromSet(self.sites[index]).write_bedgraph_file(os.path.join(self.outputDir, "bed", "sites_" + index + "_plus.bed"),  strand="+", track_options="")
            self.getGenomicArrayFromSet(self.sites[index]).write_bedgraph_file(os.path.join(self.outputDir, "bed", "sites_" + index + "_minus.bed"), strand="-", track_options="")
      else:
            self.getGenomicArrayFromSet(self.sites[index]).write_bedgraph_file(os.path.join(self.outputDir, "bed", "sites_" + index + ".bed"), strand=".", track_options="")


   #-------------------------------------------------------------
   def getGenomicArrayFromSet(self, gas):
       ga = HTSeq.GenomicArray("auto", typecode='i', stranded=self.stranded)
       for interval, values in gas.steps():
           if len(values) > 0:
               ga[interval] = len(values)

       return ga


   #-------------------------------------------------------------
   def plotOutput(self):
      outputDir = self.outputDir
      data = self.data

      print "# Plotting table for Exon Distance"


      # Plot for Exon Distances
      key = 'exondist'
      df = data[key]
      plt.plot(df.sum(axis=1))
      #plt.axis([0, 6, 0, 20])
      plt.grid(True)
      plt.xlabel('Read Length')
      plt.savefig(outputDir + '/' + key + '.png', bbox_inches='tight')

      # Plot for Intron Distances
      key = 'introndist'
      df = data[key]
      plt.plot(df.sum(axis=1))
      #ax.legend(handles, labels)
      plt.grid(True)
      plt.xlabel('Read Length')
      plt.savefig(outputDir + '/' + key + '.png', bbox_inches='tight')

      # Plot for Adapters
      #key = 'adapters'
      #df = data[key]
      #plt.plot(df.sum(axis=0))
      #plt.grid(True)
      #plt.xlabel('Adapters')
      #plt.savefig(outputDir + '/' + key + '.png', bbox_inches='tight')


      #print "# Plotting table for Intron Distance"
      #plt.plot( np.arange( -halfwinwidth, halfwinwidth ), profile )
      #plt.savefig(outputDir + '/foo.png', bbox_inches='tight')

   #-------------------------------------------------------------
   def createDir(self, directory):
       if not os.path.exists(directory):
          os.makedirs(directory)

   #-------------------------------------------------------------
   # load gtf or gff genome file of reference organism (Homo sapiens, Mus musculus) and search for exons #################
   # NOTE: for introns the gtf file has to be modified first
   def processGTF(self, gtfname):
      print "Processing the GTF file %s" % gtfname

      gtffile = HTSeq.GFF_Reader( gtfname )

      for feature in gtffile:
         if feature.type == "exon" or feature.type == "exonic":
            self.features['exons'][ feature.iv ] += feature
         if feature.type == "intron" or feature.type == "intronic":
            self.features['introns'][ feature.iv ] += feature

   #-------------------------------------------------------------
   def calcBedToFeature(self, bedname, featureUpstream='exons', featureDownstream='exons', window=399):
        print "Bed"
        bedreader = HTSeq.BED_Reader(bedname)

        for bed in bedreader:
            exon3ps = set( [ f.iv.end_d for f in self.features[featureUpstream][ bed.iv.start_d_as_pos ] ] )
            exon5ps = set( [ f.iv.start_d for f in self.features[featureDownstream][ bed.iv.start_d_as_pos ] ] )

            for e3p in exon3ps:
                dist = abs( e3p - bed.iv.end_d ) + abs(bed.iv.end_d - bed.iv.start_d)/2
                if len( exon3ps ) == 1:
                    self.increaseDataCount('bed-' + featureUpstream + '-' + featureDownstream, str( -min( dist, window ) ), "bed")


            for e5p in exon5ps:
                dist = abs( e5p - bed.iv.start_d ) + abs(bed.iv.end - bed.iv.start)/2
                if len( exon5ps ) == 1:
                    self.increaseDataCount('bed-' + featureUpstream + '-' + featureDownstream, str( min( dist, window ) ), "bed")

   #-------------------------------------------------------------

   def readExonGTF(self, gtfname):
       # store the file name
       print "Reading GTF file %s " % gtfname

       # open a gff reader connection with htseq
       gff  = HTSeq.GFF_Reader(gtfname)

       # iterate through the features of the gff file and extract the exons corresponding
       # to one sigle gene
       for geneID, exonList in self.bundleByGene( self.getExonsOnly( gff ) ):
            geneGenomicArray = HTSeq.GenomicArray( "auto", typecode="b", stranded=self.stranded )
            geneStart = 1e30
            geneEnd   = 0

            #print "-------------------------"
            # dict stores all exon ends to the start sites to filter out exons with
            # multiple start and end sites
            exonDict = {}

            # for each exon of a gene
            for exon in exonList:
                geneGenomicArray[exon.iv] = True
                geneStart    = min( geneStart, exon.iv.start )
                geneEnd      = max( geneEnd,   exon.iv.end )

                #exonDict[exon.iv.chrom exon.iv.start] = exon.iv.end

            #print "GeneStart %s" % geneStart
            #print "GeneEnd  %s " % geneEnd


            geneInterval = HTSeq.GenomicInterval( exon.iv.chrom, geneStart, geneEnd, exon.iv.strand)

            # save the gene start and the gene end
            self.increaseSiteCount('gene-start', geneInterval.start_d_as_pos, geneID)
            self.increaseSiteCount('gene-end',   geneInterval.end_d_as_pos, geneID)

            #print "[ %s ] -------  " % geneID
            mem = True

            #print geneGenomicArray
            i = 1
            e = 1
            for interval, isExon in geneGenomicArray[ geneInterval ].steps():

                # print  '%s -- %s ' % (interval, isExon)
                if isExon:
                    #add = True
                    #if self.filterExons:
                        #if
                    #TODO Filter Exons

                    # print "exon %s " % str(interval)
                    self.features['exons'][interval] = str(geneID) + "-e" + str(e)
                    e += 1

                    ## add intron exon site
                    #if(geneStart == interval.start):
                       # print "first Exon"

                    #if(geneEnd == interval.end):
                       # print "last Exon"

                    if mem:
                        mem = False
                    else:
                        sys.os._exit(1)

                else:
                    # Stores the correct exon intron site (site on the intron)
                    startPos = 0
                    endPos = 0

                    self.increaseSiteCount('exon-intron', interval.start_d_as_pos, geneID)
                    # Sotres the correct specific intron exon site
                    self.increaseSiteCount('intron-exon', interval.end_d_as_pos, geneID)
                    #print "exon-intron"
                    #print interval.start_d_as_pos
                    #print "intron-exon"
                    #print interval.end_d_as_pos
                    #print "intron %s " % str(interval)
                    #print "exon-intron %s " % str(interval.start_d_as_pos)
                    #print "intron-exon %s " % str(interval.end_d_as_pos)
                    self.features['introns'][interval] = str(geneID) + "-i" + str(i)
                    i += 1

                    if not mem:
                        mem = True
                    else:
                        sys.os._exit(1)


   #-------------------------------------------------------------
   def getTimeStamp(self):
       return time.strftime("%c")

   #-------------------------------------------------------------
   # returns exons Only
   def getExonsOnly( self, features ):
        for feature in features:
            if feature.type == "exon":
                yield feature

   #-------------------------------------------------------------
   # Bundles the features by gene
   def bundleByGene( self, features ):
        currentID = None
        genelist = []
        for feature in features:
            if feature.attr["gene_id"] == currentID:
                genelist.append( feature )
            else:
                if currentID is not None:
                    yield ( currentID, genelist )
                genelist = [ feature ]
                currentID = feature.attr["gene_id"]

        yield(currentID, genelist)


