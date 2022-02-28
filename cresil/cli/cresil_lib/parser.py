import re, os
        	
class refgene_render:
    """handle bed format
       input in list or str
    """
    def __init__(self, refgene_in, from_string=False, bins=True):
        refgene = []
        if from_string or type(refgene_in) == str:
            refgene = re.split("\s+", refgene_in.strip())
        else:
            refgene = refgene_in
        if bins: refgene = refgene[1:]
        if not len(refgene) ==  15:
            print ("wrong refgene format (len={}) of '{}'".format(len(refgene),refgene_in))
        else:            
            self.name = refgene[0]
            self.chrom = refgene[1]
            self.strand = refgene[2]
            self.txStart = int(refgene[3])
            self.txEnd = int(refgene[4])
            self.cdsStart = int(refgene[5])
            self.cdsEnd = int(refgene[6])
            self.exonCount = int(refgene[7])
            self.exonStarts = refgene[8]
            self.exonEnds = refgene[9]
            self.score = int(refgene[10])
            self.symbol = refgene[11]
            self.cdsStartStat = refgene[12]
            self.cdsEndStat = refgene[13]
            self.exonFrames = refgene[14]
            self.len = len(refgene)
            
    def __len__(self):
        return self.len

    def get_exon_starts(self):
        exonStarts = map(int,self.exonStarts.strip(",").split(","))
        return exonStarts
    
    def get_exon_ends(self):
        exonEnds = map(int,self.exonEnds.strip(",").split(","))
        return exonEnds
    
    def get_exon_frames(self):
        exonFrames = map(int,self.exonFrames.strip(",").split(","))
        return exonFrames

    def get_bed(self, bed_format=12, symbol_name=False, as_string=False):
        chrom = self.chrom
        start = self.txStart
        end = self.txEnd
        if symbol_name:
            name = self.symbol
        else:
            name = self.name
        score = self.score
        strand = self.strand
        color = "0"
        thickStart = self.cdsStart
        thickEnd = self.cdsEnd
        blocks = self.exonCount
        blockSizes = [ e-s for s, e in zip(self.get_exon_starts(), self.get_exon_ends())]
        blockStarts = [ x - self.txStart for x in self.get_exon_starts()]

        if as_string:
            return "\t".join(map(str,[chrom, start, end, name, score, strand, thickStart, thickEnd, color, blocks, "{},".format(",".join(map(str,blockSizes))), "{},".format(",".join(map(str,blockStarts)))]))
        else:
            return [chrom, start, end, name, score, strand, thickStart, thickEnd, color, blocks, "{},".format(",".join(map(str,blockSizes))), "{},".format(",".join(map(str,blockStarts)))]

    def get_gff(self, as_string=False, as_list=True):
        out_location = []
        exonFrames = self.get_exon_frames()
        exonStarts = self.get_exon_starts()
        exonEnds = self.get_exon_ends()

        ## To fix muti-start and stop codon
        for i in range(self.exonCount):
            out_location.append(['exon', exonStarts[i]+1, exonEnds[i], '.'])
            if exonFrames[i] >= 0:
                if self.cdsStart > exonStarts[i] and self.cdsStart < exonEnds[i]:
                    if self.strand == "-":
#                        out_location.append(['stop_codon', self.cdsStart+1, self.cdsStart+3, '.'])
                        out_location.append(['CDS', self.cdsStart+1, exonEnds[i], exonFrames[i]])
                    else:
#                        out_location.append(['start_codon', self.cdsStart+1, self.cdsStart+3, '.'])
                        out_location.append(['CDS', self.cdsStart+1, exonEnds[i], exonFrames[i]])
                elif self.cdsEnd > exonStarts[i] and self.cdsEnd < exonEnds[i]:
                    if self.strand == "-":
                        out_location.append(['CDS', exonStarts[i]+1, self.cdsEnd, exonFrames[i]])
#                        out_location.append(['start_codon', self.cdsEnd-2, self.cdsEnd, '.'])
                    else:
                        out_location.append(['CDS', exonStarts[i]+1, self.cdsEnd, exonFrames[i]])
#                        out_location.append(['stop_codon', self.cdsEnd-3+1, self.cdsEnd, '.'])
                else:
                    out_location.append(['CDS', exonStarts[i]+1, exonEnds[i], exonFrames[i]])
        
        gff = []
        for anno, start, end, frame in out_location:
##            gff.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id "%s"; gene_name="%s"; transcript_id="%s";'\
##                       %(self.chrom, "refseq_render", anno, start, end, ".", self.strand,\
##                         frame, self.symbol, self.symbol, self.name))
            gff.append([self.chrom, "refseq_render", anno, start, end, ".",self.strand, frame,\
                        'gene_id "{}"; transcript_id "{}";'.format(self.symbol, self.name)])
        if as_string:
            return "\n".join(["\t".join(map(str,L)) for L in gff])
        elif as_list:
            return ["\t".join(map(str,L)) for L in gff]
        else:
            return gff
        

class refFlat(refgene_render):
    def __init__(self, refgene_in, from_string=False):
        refgene = []
        if from_string:
            refgene = re.split("\s+", refgene_in.strip())
        else:
            refgene = refgene_in

        if not len(refgene) ==  11:
            print (f"wrong refFlat format of '{refgene_in}'" )
        else:            
            self.symbol = refgene[0]
            self.name = refgene[1]
            self.chrom = refgene[2]
            self.strand = refgene[3]
            self.txStart = int(refgene[4])
            self.txEnd = int(refgene[5])
            self.cdsStart = int(refgene[6])
            self.cdsEnd = int(refgene[7])
            self.exonCount = int(refgene[8])
            self.exonStarts = refgene[9]
            self.exonEnds = refgene[10]

            self.score = "0"
            self.cdsStartStat = False
            self.cdsEndStat = False
            self.exonFrames = False

            self.len = len(refgene)



class bed_render:
    """handle bed format
       input in list or str
    """
    import re
    def __init__(self, bed_in="chr1	901876	910484	NM_032129	0	+	908364	909361	0	16	118,100,147,81,73,128,96,81,76,137,150,141,141,219,49,663,	0,207,3780,4024,4189,4382,4616,4827,5578,5791,6364,6689,7003,7336,7819,7945,", bedtype=12, from_string=True):
        bed = []
        if from_string and type(bed_in) == str:
            bed = re.split("\s+", bed_in.strip())
        else:
            bed = bed_in
        
        if len(bed) < bedtype: 
            if len(bed)>=6:
                bedtype = 6
            elif len(bed)<6 and len(bed)>=4:
                bedtype = 4
            elif len(bed)<6 and len(bed)>=3:
                bedtype = 3
            else:
                print ("Ambiguous BED and BED type")
        self.bedtype = bedtype
        self.bed = bed[:bedtype]
        self.bedExtra = bed[bedtype:]
        self.strand = "+"
        self.name = "."
        self.score = "0"
        if len(bed) >= 3:
            self.chrom = bed[0]
            self.start = int(bed[1])
            self.end = int(bed[2])
            self.bed[1] = self.start
            self.bed[2] = self.end
            self.len = self.end-self.start
        if len(bed) >= 4:
            self.name = bed[3]
        if len(bed) >= 6:
            self.score = bed[4]
            self.bed[4] = self.score
            self.strand = bed[5]
            self.tss = 0
            if self.strand == "-":
                self.tss = self.end
                self.tes = self.start
            else:
                self.tss = self.start
                self.tes = self.end
        if self.strand != "+" and self.strand != "-":
            self.strand = "+"
                
        if len(bed) >= 12:
            self.thickStart = int(bed[6])
            self.thickEnd = int(bed[7])
            self.bed[6] = self.thickStart
            self.bed[7] = self.thickEnd
            self.color = bed[8]
            self.blocks = int(bed[9])
            self.bed[9] = self.blocks
            self.blockSizes = bed[10]
            self.blockStarts = bed[11]


    def get_chrom_starts(self):
        blockStarts = [int(x) for x in self.blockStarts.strip(",").split(",")]
        return [x+self.start for x in blockStarts]
    
    def get_chrom_ends(self):
        blockSizes = [int(x) for x in self.blockSizes.strip(",").split(",")]
        Starts = self.get_chrom_starts()
        return [Starts[i]+blockSizes[i] for i in range(self.blocks)]

    def get_location_specific_block(self, block, as_string=False, separator=" "):
        if block > self.blocks or block < 1:
            print ("block out of range")
            return False
        starts = self.get_chrom_starts()
        ends = self.get_chrom_ends()
        if as_string:
            return separator.join(map(str,[self.chrom, starts[block-1], ends[block-1]]))
        else:
            return [self.chrom, starts[block-1], ends[block-1]]

    def get_gff(self, as_string=False, as_list=True):
        out_location = []
        exonStarts = self.get_chrom_starts()
        exonEnds = self.get_chrom_ends()

        ## To fix muti-start and stop codon
        for i in range(self.blocks):
            out_location.append(['exon', exonStarts[i]+1, exonEnds[i], '.'])
            
            if self.thickStart != self.thickEnd and exonStarts[i] < self.thickEnd and exonEnds[i] > self.thickStart:
                if self.thickStart > exonStarts[i] and self.thickStart < exonEnds[i]:
                    out_location.append(['CDS', self.thickStart+1, exonEnds[i], '.'])
                elif self.thickEnd > exonStarts[i] and self.thickEnd < exonEnds[i]:
                    out_location.append(['CDS', exonStarts[i]+1, self.thickEnd, '.'])
                else:
                    out_location.append(['CDS', exonStarts[i]+1, exonEnds[i], '.'])
        
        gff = []
        for anno, start, end, frame in out_location:
##            gff.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id "%s"; gene_name="%s"; transcript_id="%s";'\
##                       %(self.chrom, "refseq_render", anno, start, end, ".", self.strand,\
##                         frame, self.symbol, self.symbol, self.name))
            gff.append([self.chrom, "bed_render", anno, start, end, ".",self.strand, frame,\
                        'gene_id "{}"; transcript_id "{}";'.format(self.name, self.name)])
        if as_string:
            return "\n".join(["\t".join(map(str,L)) for L in gff])
        elif as_list:
            return ["\t".join(map(str,L)) for L in gff]
        else:
            return gff

    def get_ucsc_refseq(self, gene="", as_string=False, as_list=True):
        bin_t = 0
        name = self.name
        if gene == "": gene = name
        chrom = self.chrom
        strand = self.strand
        start = self.start
        end = self.end
        cdsStart = self.thickStart
        cdsEnd = self.thickEnd
        exonCount = self.blocks
        exonStarts = "{},".format(",".join(map(str,self.get_chrom_starts())))
        exonEnds = "{},".format(",".join(map(str,self.get_chrom_ends())))
        score = self.score
        name2 = gene
        cdsStartStat = "unk"
        cdsEndStat = "unk"
        exonFrames = "{},".format(",".join(['0']*exonCount))

        if as_string:
            return "\t".join(map(str,[bin_t,name,chrom,strand,start,end,cdsStart,cdsEnd,exonCount,exonStarts, \
                                             exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames]))
        else:
            return [bin_t,name,chrom,strand,start,end,cdsStart,cdsEnd,exonCount,exonStarts, \
                    exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames]

    def get_refFlat(self, gene="", as_string=False, as_list=True):
        name = self.name
        if gene == "": gene = name
        chrom = self.chrom
        strand = self.strand
        start = self.start
        end = self.end
        cdsStart = self.thickStart
        cdsEnd = self.thickEnd
        exonCount = self.blocks
        exonStarts = "{},".format(",".join(map(str,self.get_chrom_starts())))
        exonEnds = "{},".format(",".join(map(str,self.get_chrom_ends())))
        name2 = gene

        if as_string:
            return "\t".join(map(str,[name2,name,chrom,strand,start,end,cdsStart,cdsEnd,exonCount,exonStarts, \
                                             exonEnds]))
        else:
            return [name2,name,chrom,strand,start,end,cdsStart,cdsEnd,exonCount,exonStarts, \
                    exonEnds]

    def stitch_bed(self, bed, as_string=False):
        chrom = self.chrom
        start = bed[0][1]
        end = bed[-1][2]
##        name = ".".join(bed[0][3].split(".")[:2])
        name = bed[0][3]
        score = self.score
        strand = self.strand
        color = "0"
        thickStart = start
        thickEnd = end
        blocks = len(bed)
        blockSizes = []
        blockStarts = []
        for L in bed:
            s = int(L[1])
            e = int(L[2])
            blockSizes.append( e-s )
            blockStarts.append( s - start )
        if as_string:
            return "{}\n".format("\t".join(map(str,[chrom, start, end, name, score, strand, thickStart, thickEnd, color, blocks, "{},".format(",".join(map(str,blockSizes))), "{},".format(",".join(map(str,blockStarts)))])))
        else:
            return [chrom, start, end, name, score, strand, thickStart, thickEnd, color, blocks, "{},".format(",".join(map(str,blockSizes))), "{},".format(",".join(map(str,blockStarts)))]
        

    def get_feature(self, feature="exon", as_string=False, stitch=False, verbose=False):
        """  features = [exon, intron, cds, 5utr, 3utr] """
##        print "Hello"
        bed = []
        if self.bedtype == 12:
            if feature == "exon":
                starts = self.get_chrom_starts()
                ends = self.get_chrom_ends()
                exon_first_last = ""
                if self.strand == '-':
                    if self.blocks == 1:
                        exon_first_last = ".single"
                        if starts[0] > ends[0]:
                            if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[0],starts[0]))
                            return False
                        if stitch:
                            bed.append([self.chrom, starts[0], ends[0], "{}.exon".format(self.name), self.score, self.strand])
                        else:
                            bed.append([self.chrom, starts[0], ends[0], "{}.exon.{}{}".format(self.name,int(self.blocks)-0,exon_first_last), self.score, self.strand])
                    else:
                        for i in range(self.blocks):
                            if i == 0:
                                exon_first_last = ".l"
                            elif i+1 == int(self.blocks):
                                exon_first_last = ".f"
                            else:
                                exon_first_last = ""
                            if starts[i] > ends[i]:
                                if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[0],starts[0]))
                                return False
                            if stitch:
                                bed.append([self.chrom, starts[i], ends[i], "{}.exon".format(self.name), self.score, self.strand])
                            else:
                                bed.append([self.chrom, starts[i], ends[i], "{}.exon.{}{}".format(self.name,int(self.blocks)-i,exon_first_last), self.score, self.strand])
                else:
                    for i in range(self.blocks):
                        if i == 0:
                            exon_first_last = ".f"
                        elif i+1 == int(self.blocks):
                            exon_first_last = ".l"
                        else:
                            exon_first_last = ""
                        if starts[i] > ends[i]:
                            if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[i],starts[i]))
                            return False
                        if stitch:
                            bed.append([self.chrom, starts[i], ends[i], "{}.exon".format(self.name), self.score, self.strand])
                        else:
                            bed.append([self.chrom, starts[i], ends[i], "{}.exon.{}{}".format(self.name,i+1,exon_first_last), self.score, self.strand])
                if as_string:
                    if stitch:
                        return self.stitch_bed(bed, as_string=True)
                    else:                    
                        return "{}\n".format("".join(["\t".join(map(str,x))+"\n" for x in bed]))
                else:
                    if stitch:
                        return self.stitch_bed(bed)
                    else:
                        return bed
            
            elif feature == "intron":
                bed = []
                if self.blocks > 1:
                    starts = self.get_chrom_ends()[:-1]
                    ends = self.get_chrom_starts()[1:]
                    first_last = ""
                    if self.strand == '-':
                        if len(starts) == 1:
                            first_last = ".single"
                            if starts[0] > ends[0]:
                                if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[0],starts[0]))
                                return False
                            if stitch:
                                bed.append([self.chrom, starts[0], ends[0], "{}.intron".format(self.name), self.score, self.strand])
                            else:
                                bed.append([self.chrom, starts[0], ends[0], "{}.intron.{}{}".format(self.name,len(starts)-0,first_last), self.score, self.strand])
                        else:
                            for i in range(len(starts)):
                                if i == 0:
                                    first_last = ".l"
                                elif i+1 == int(len(starts)):
                                    first_last = ".f"
                                else:
                                    first_last = ""
                                if starts[i] > ends[i]:
                                    if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[i],starts[i]))
                                    return False
                                if stitch:
                                    bed.append([self.chrom, starts[i], ends[i], "{}.intron".format(self.name), self.score, self.strand])
                                else:
                                    bed.append([self.chrom, starts[i], ends[i], "{}.intron.{}{}".format(self.name,len(starts)-i,first_last), self.score, self.strand])
                    else:
                        for i in range(len(starts)):
                            if i == 0:
                                first_last = ".f"
                            elif i+1 == int(len(starts)):
                                first_last = ".l"
                            else:
                                first_last = ""
                            if starts[i] > ends[i]:
                                if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[i],starts[i]))
                                return False
                            if stitch:
                                bed.append([self.chrom, starts[i], ends[i], "{}.intron".format(self.name), self.score, self.strand])
                            else:
                                bed.append([self.chrom, starts[i], ends[i], "{}.intron.{}{}".format(self.name,i+1,first_last), self.score, self.strand])
                else:
                    if verbose: 
                        print (f"[Error] no Intron found in {self.name}")
                    return False
                if as_string:
                    if stitch:
                        return self.stitch_bed(bed, as_string=True)
                    else:                    
                        return "{}\n".format("".join(["\t".join(map(str,x))+"\n" for x in bed]))
                else:
                    if stitch:
                        return self.stitch_bed(bed)
                    else:
                        return bed

            elif feature == "cds":
##                self.thickStart
##                self.thickEnd
                if self.thickStart == self.thickEnd:
                    if verbose: print ("no CDS regions")
                    return False
                starts = list(filter(lambda x: x > self.thickStart and x < self.thickEnd, self.get_chrom_starts()))
                starts.insert(0, self.thickStart)
                ends = list(filter(lambda x: x < self.thickEnd and x > self.thickStart, self.get_chrom_ends()))
                ends.append(self.thickEnd)
##                print len(starts), starts
##                print len(ends), ends
                exon_first_last = ""
                if len(starts) != len(ends):
                    if verbose: print ("[Error] number of starts and ends are not equal")
                    return False
##                    break
                if self.strand == '-':
                    if len(starts) == 1:
                        exon_first_last = ".single"
                        if starts[0] > ends[0]:
                            if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[0],starts[0]))
                            return False
                        if stitch:
                            bed.append([self.chrom, starts[0], ends[0], "{}.cds".format(self.name), self.score, self.strand])
                        else:
                            bed.append([self.chrom, starts[0], ends[0], "{}.cds.{}{}".format(self.name,int(len(starts))-0,exon_first_last), self.score, self.strand])
                    else:
                        for i in range(len(starts)):
                            if i == 0:
                                exon_first_last = ".l"
                            elif i+1 == int(len(starts)):
                                exon_first_last = ".f"
                            else:
                                exon_first_last = ""
                            if starts[i] > ends[i]:
                                if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[i],starts[i]))
                                return False
                            if stitch:
                                bed.append([self.chrom, starts[i], ends[i], "{}.cds".format(self.name), self.score, self.strand])
                            else:
                                bed.append([self.chrom, starts[i], ends[i], "{}.cds.{}{}".format(self.name,int(len(starts))-i,exon_first_last), self.score, self.strand])
                else:
                    for i in range(len(starts)):
                        if i == 0:
                            exon_first_last = ".f"
                        elif i+1 == int(len(starts)):
                            exon_first_last = ".l"
                        else:
                            exon_first_last = ""
                        if starts[i] > ends[i]:
                                if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[i],starts[i]))
                                return False
                        if stitch:
                            bed.append([self.chrom, starts[i], ends[i], "{}.cds".format(self.name), self.score, self.strand])
                        else:
                            bed.append([self.chrom, starts[i], ends[i], "{}.cds.{}{}".format(self.name,i+1,exon_first_last), self.score, self.strand])
                if as_string:
                    if stitch:
                        return self.stitch_bed(bed, as_string=True)
                    else:                    
                        return "{}\n".format("".join(["\t".join(map(str,x))+"\n" for x in bed]))
                else:
                    if stitch:
                        return self.stitch_bed(bed)
                    else:
                        return bed

            elif feature in ['5utr','3utr']:
##                self.thickStart
##                self.thickEnd
                if self.thickStart == self.thickEnd:
                    if verbose: print ("[Error] no CDS regions")
                    return False
                
                # if self.thickStart >= x --> fail
                L_starts = []
                for x in self.get_chrom_starts():
                    if x < self.thickStart:
                        L_starts.append(x)
                L_ends  = []
                for x in self.get_chrom_ends():
                    if x < self.thickStart:
                        L_ends.append(x)
                L_ends.append(self.thickStart)
         
                R_starts = [self.thickEnd]
                for x in self.get_chrom_starts():
                    # print(f"{x} > {self.thickEnd} = {x > self.thickEnd}")
                    if x > self.thickEnd:
                        R_starts.append(x)

                R_ends = []
                for x in self.get_chrom_ends():
                    # print(f"{x} > {self.thickEnd} = {x > self.thickEnd}")
                    if x > self.thickEnd:
                        R_ends.append(x)

                if self.strand == '-':
                    if feature == '3utr':
                        starts = L_starts
                        ends = L_ends
                    else:
                        starts = R_starts
                        ends = R_ends
                else:
                    if feature == '5utr':
                        starts = L_starts
                        ends = L_ends
                    else:
                        starts = R_starts
                        ends = R_ends

                if len(starts) == 0 or len(ends)==0:
                    if verbose: print (f"no {feature}" )
                    return False
                if len(starts) != len(ends):
                    if verbose: print ("[Error] start points and end points are not equal")
                    return False


                for i in range(len(starts)):
                    if starts[i] > ends[i]:
                        if verbose: print ("[Error] end ({}) need to be higher than start ({})".format(ends[i],starts[i]))
                        return False
                    bed.append([self.chrom, starts[i], ends[i], "{}.{}".format(self.name,feature), self.score, self.strand])

                if as_string:
                    if stitch:
                        return self.stitch_bed(bed, as_string=True)
                    else:                    
                        return "{}\n".format("".join(["\t".join(map(str,x))+"\n" for x in bed]))
                else:
                    if stitch:
                        return self.stitch_bed(bed)
                    else:
                        return bed

        else:
            if verbose: print ("this module require BED12")
            return False
##            if as_string:
##                if stitch:
##                    return self.stitch_bed([self.bed], as_string=True)
##                else:                    
##                    return "%s\n"%"".join(["\t".join(map(str,x))+"\n" for x in [self.bed]])
##            else:
##                if stitch:
##                    return self.stitch_bed([self.bed])
##                else:
##                    return self.bed
##
##            if as_string:
##                return "%s\n"%"".join(["\t".join(map(str,x))+"\n" for x in [self.bed]])
##            else:
##                return self.bed
    
    def get_bed(self, format_num=12, as_string=False, separator="\t"):
        bed = self.bed
        if format_num == 12 and len(bed) == 12:
            if as_string:
                return  ("{}\n".format(separator.join(map(str,bed))))
            else:
                return (bed)
        elif format_num == 6 or len(bed) == 6:
            if as_string:
                return  ("{}\n".format(separator.join(map(str,bed[:6]))))
            else:
                return (bed[:6])
        elif format_num == 4 or len(bed) == 4:
            if as_string:
                return  ("{}\n".format(separator.join(map(str,bed[:4]))))
            else:
                return( bed[:3])
        elif format_num == 3 or len(bed) == 3:
            if as_string:
                return(  "{}\n".format(separator.join(map(str,bed[:3]))))
            else:
                return( bed[:3])
        else:
            print ("Please use only 3, 4, 6, and 12")
            
    def get_flanks(self, up=0, down=0, region="gene", genome="", not_fit_genome_size=False, stranded=True, as_string=False):
        ## add genome location file.
#        if len(chromosome_sizes)==0:
#            print "does not have chromosome_sizes information"
 #           return dict()
        if genome != "":
            if os.path.isfile(genome):           
                chromosome_sizes = {}
                for l in open(genome).xreadlines():
                    L = l.strip().split("\t")
                    chromosome_sizes[L[0]] = (0,int(L[1]))
            else:
                exit(f'There is no genome assembly or "{genome}" is not exist')

        # flanks = {}
        if stranded and self.strand == '-':
            if region == "gene":
                start = int(self.start)-down
                end = int(self.end)+up
            elif region == "tss":
                start = int(self.end)-down
                end = int(self.end)+up
            elif region == "tes":
                start = int(self.start)-down
                end = int(self.start)+up
        elif stranded and self.strand == '+':
            if region == "gene":
                start = int(self.start)-up
                end = int(self.end)+down
            elif region == "tss":
                start = int(self.start)-up
                end = int(self.start)+down
            elif region == "tes":
                start = int(self.end)-up
                end = int(self.end)+down
        else:
            if region == "gene":
                start = int(self.start)-up
                end = int(self.end)+down
            elif region == "tss":
                start = int(self.start)-up
                end = int(self.start)+down
            elif region == "tes":
                start = int(self.end)-up
                end = int(self.end)+down
            
        name = self.name
        up_name = up
        down_name = down
        if genome != "" : 
            fit_genome = True        
            if start < chromosome_sizes[self.chrom][0]:
                start = chromosome_sizes[self.chrom][0]
                up = self.start - chromosome_sizes[self.chrom][0]
                fit_genome = False
            if end > chromosome_sizes[self.chrom][1]:
                end = chromosome_sizes[self.chrom][1]
                down = chromosome_sizes[self.chrom][1] - self.end
                fit_genome = False
        else:
            fit_genome = False
        name = "{}_{}_{}up_{}down".format(name,region,up_name,down_name)

#        if start < 0 or end > chromosome_sizes[self.chrom]:
#            flanks['start'] = start
#            flanks['end'] = end
#            flanks['size'] = -1
#        else:
#            flanks['start'] = start
#            flanks['end'] = end
#            flanks['size'] = end - start
        if fit_genome:
            if as_string:
                return "{}\n".format("\t".join(map(str, [self.chrom, start, end, name, self.score, self.strand])))
            else:
                return [self.chrom, start, end, name, self.score, self.strand]
        elif not_fit_genome_size:
            if as_string:
                return "{}\n".format("\t".join(map(str, [self.chrom, start, end, name, self.score, self.strand])))
            else:
                return [self.chrom, start, end, name, self.score, self.strand]
        else:
            return False
    
    def set_promoter(self, up=2000, down=2000, region="tss"):
        self.promoter = self.get_flanks(up, down, region="tss")

    def set_downstream(self, up=2000, down=2000, region="tss"):
        self.downstream = self.get_flanks(up, down, region="tes")

                
