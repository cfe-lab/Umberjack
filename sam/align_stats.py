class AlignStats:
    def __init__(self):
        self.total_insert_1mate_hi_qual = 0  # high quality inserted base in nonoverlapping position of mate
        self.total_insert_1mate_lo_qual = 0  # low quality inserted base in nonoverlapping position of mate
        self.total_insert_1mate = 0  # inserted base in nonoverlapping position of mate
        self.total_insert_agree = 0 # inserted base in overlapping position of mates, base shared by both mates
        self.total_insert_agree_hi_qual = 0  # high quality inserted base in overlapping position of mates, base shared by both mates
        self.total_insert_agree_lo_qual = 0  # low quality inserted base in overlapping position of mates, base shared by both mates
        self.total_insert_conflict = 0  # inserted base in overlapping position of mates, base not shared by both mates
        self.total_insert_conflict_hi_qual = 0  # high quality inserted base in overlapping position of mates, base not shared by both mates
        self.total_insert_conflict_hilo_qual = 0  # inserted base in overlapping position of mates, base not shared by both mates.  1 high qual, 1 low qual
        self.total_insert_conflict_lo_qual = 0  # low quality inserted base in overlapping position of mates, base not shared by both mates
        self.total_insert_conflict_hino_qual = 0  # inserted base in overlapping position of mates.  1 high quality, 1 not existing.
        self.total_insert_conflict_lono_qual = 0  # inserted base in overlapping position of mates.  1 low quality, 1 not existing.
        self.total_insert_conflict_equal_hi_qual = 0  # inserted base in overlapping position of mates.  base not shared by both mates. both same high quality.
        self.total_inserts = 0  # total inserted bases
        self.total_insert_blocks = 0  #  total blocks of contiguous inserted bases

        # TODO:  keep histo of insert length size

        self.total_match_1mate = 0  #  base in nonoverlapping position of mate
        self.total_match_1mate_hi_qual = 0  #  high quality base in nonoverlapping position of mate
        self.total_match_1mate_lo_qual = 0  #  low quality base in nonoverlapping position of mate
        self.total_match_nonconflict = 0  # base in overlapping position of mates, base shared by both mates
        self.total_match_nonconflict_hi_qual = 0  # high quality base in overlapping position of mates, base shared by both mates
        self.total_match_nonconflict_lo_qual = 0  # low quality base in overlapping position of mates, base shared by both mates
        self.total_match_conflict = 0  # base in overlapping position of mates, base not shared by both mates
        self.total_match_conflict_lo_qual = 0  # bases in overlapping position of mates, bases not shared by both mates.  Both mates low qual.
        self.total_match_conflict_hilo_qual = 0  # base not shared by both mates, 1 high quality, 1 low quality
        self.total_match_conflict_hi_qual = 0  # bases in overlapping position of mates, bases not shared by both mates.  both mates high quality
        self.total_match_conflict_equal_hi_qual = 0  # bases in overlapping position of mates, bases not shared by both mates.  both mates high qual



    def merge_stats(self, other_stats):
        self.total_insert_1mate_hi_qual += other_stats.total_insert_1mate_hi_qual
        self.total_insert_1mate_lo_qual += other_stats.total_insert_1mate_lo_qual
        self.total_insert_1mate += other_stats.total_insert_1mate
        self.total_insert_agree += other_stats.total_insert_agree
        self.total_insert_agree_hi_qual += other_stats.total_insert_agree_hi_qual # high quality inserted base in overlapping position of mates, base shared by both mates
        self.total_insert_agree_lo_qual += other_stats.total_insert_agree_lo_qual  # low quality inserted base in overlapping position of mates, base shared by both mates
        self.total_insert_conflict += other_stats.total_insert_conflict  # inserted base in overlapping position of mates, base not shared by both mates
        self.total_insert_conflict_hi_qual += other_stats.total_insert_conflict_hi_qual  # high quality inserted base in overlapping position of mates, base not shared by both mates
        self.total_insert_conflict_hilo_qual += other_stats.total_insert_conflict_hilo_qual  # inserted base in overlapping position of mates, base not shared by both mates.  1 high qual, 1 low qual
        self.total_insert_conflict_lo_qual += other_stats.total_insert_conflict_lo_qual  # low quality inserted base in overlapping position of mates, base not shared by both mates
        self.total_insert_conflict_hino_qual += other_stats.total_insert_conflict_hino_qual  # inserted base in overlapping position of mates.  1 high quality, 1 not existing.
        self.total_insert_conflict_lono_qual += other_stats.total_insert_conflict_lono_qual  # inserted base in overlapping position of mates.  1 low quality, 1 not existing.
        self.total_insert_conflict_equal_hi_qual += other_stats.total_insert_conflict_equal_hi_qual  # bases in overlapping position of mates, bases not shared by both mates.  both mates high qual
        self.total_inserts += other_stats.total_inserts  # total inserted bases
        self.total_insert_blocks += other_stats.total_insert_blocks  #  total blocks of contiguous inserted bases


        self.total_match_1mate += other_stats.total_match_1mate  #  base in nonoverlapping position of mate
        self.total_match_1mate_hi_qual += other_stats.total_match_1mate_hi_qual  #  high quality base in nonoverlapping position of mate
        self.total_match_1mate_lo_qual += other_stats.total_match_1mate_lo_qual  #  low quality base in nonoverlapping position of mate
        self.total_match_nonconflict += other_stats.total_match_nonconflict  # base in overlapping position of mates, base shared by both mates
        self.total_match_nonconflict_hi_qual += other_stats.total_match_nonconflict_hi_qual  # high quality base in overlapping position of mates, base shared by both mates
        self.total_match_nonconflict_lo_qual += other_stats.total_match_nonconflict_lo_qual  # low quality base in overlapping position of mates, base shared by both mates
        self.total_match_conflict += other_stats.total_match_conflict  # base in overlapping position of mates, base not shared by both mates
        self.total_match_conflict_lo_qual += other_stats.total_match_conflict_lo_qual  # bases in overlapping position of mates, bases not shared by both mates.  Both mates low qual.
        self.total_match_conflict_hilo_qual += other_stats.total_match_conflict_hilo_qual  # bases in overlapping position of mates, bases not shared by both mates. 1 high qual, 1 low qual.
        self.total_match_conflict_hi_qual += other_stats.total_match_conflict_hi_qual  # bases in overlapping position of mates, bases not shared by both mates.  both mates high qual
        self.total_match_conflict_equal_hi_qual += other_stats.total_match_conflict_equal_hi_qual  # bases in overlapping position of mates, bases not shared by both mates.  both mates high qual

    def dump_stats(self):
        output = ("total_insert_1mate_hi_qual={}\n".format(self.total_insert_1mate_hi_qual) +
                  "total_insert_1mate_lo_qual={}\n".format(self.total_insert_1mate_lo_qual) +
                  "total_insert_1mate={}\n".format(self.total_insert_1mate) +
                  "total_insert_agree={}\n".format(self.total_insert_agree) +
                  "total_insert_agree_hi_qual={}\n".format(self.total_insert_agree_hi_qual) +
                  "total_insert_agree_lo_qual={}\n".format(self.total_insert_agree_lo_qual) +
                  "total_insert_conflict={}\n".format(self.total_insert_conflict) +
                  "total_insert_conflict_hi_qual={}\n".format(self.total_insert_conflict_hi_qual) +
                  "total_insert_conflict_hilo_qual={}\n".format(self.total_insert_conflict_hilo_qual) +
                  "total_insert_conflict_lo_qual={}\n".format(self.total_insert_conflict_lo_qual) +
                  "total_insert_conflict_hino_qual={}\n".format(self.total_insert_conflict_hino_qual) +
                  "total_insert_conflict_lono_qual={}\n".format(self.total_insert_conflict_lono_qual) +
                  "total_insert_conflict_equal_hi_qual={}\n".format(self.total_insert_conflict_equal_hi_qual) +
                  "total_inserts={}\n".format(self.total_inserts) +
                  "total_insert_blocks={}\n".format(self.total_insert_blocks) +
                  "total_match_1mate={}\n".format(self.total_match_1mate) +
                  "total_match_1mate_hi_qual={}\n".format(self.total_match_1mate_hi_qual) +
                  "total_match_1mate_lo_qual={}\n".format(self.total_match_1mate_lo_qual) +
                  "total_match_nonconflict={}\n".format(self.total_match_nonconflict) +
                  "total_match_nonconflict_hi_qual={}\n".format(self.total_match_nonconflict_hi_qual) +
                  "total_match_nonconflict_lo_qual={}\n".format(self.total_match_nonconflict_lo_qual) +
                  "total_match_conflict={}\n".format(self.total_match_conflict) +
                  "total_match_conflict_lo_qual={}\n".format(self.total_match_conflict_lo_qual) +
                  "total_match_conflict_hilo_qual={}\n".format(self.total_match_conflict_hilo_qual) +
                  "total_match_conflict_hi_qual={}\n".format(self.total_match_conflict_hi_qual) +
                  "total_match_conflict_equal_hi_qual={}\n".format(self.total_match_conflict_equal_hi_qual))

        return output


    def dump_insert_stats(self):
        output = ("total_insert_1mate_hi_qual={} ".format(self.total_insert_1mate_hi_qual) +
                  "total_insert_1mate_lo_qual={} ".format(self.total_insert_1mate_lo_qual) +
                  "total_insert_1mate={} ".format(self.total_insert_1mate) +
                  "total_insert_agree={} ".format(self.total_insert_agree) +
                  "total_insert_agree_hi_qual={} ".format(self.total_insert_agree_hi_qual) +
                  "total_insert_agree_lo_qual={} ".format(self.total_insert_agree_lo_qual) +
                  "total_insert_conflict={} ".format(self.total_insert_conflict) +
                  "total_insert_conflict_hi_qual={} ".format(self.total_insert_conflict_hi_qual) +
                  "total_insert_conflict_hilo_qual={} ".format(self.total_insert_conflict_hilo_qual) +
                  "total_insert_conflict_lo_qual={} ".format(self.total_insert_conflict_lo_qual) +
                  "total_insert_conflict_hino_qual={} ".format(self.total_insert_conflict_hino_qual) +
                  "total_insert_conflict_lono_qual={}".format(self.total_insert_conflict_lono_qual) +
                  "total_insert_conflict_equal_hi_qual={}\n".format(self.total_insert_conflict_equal_hi_qual) +
                  "total_inserts={} ".format(self.total_inserts) +
                  "total_insert_blocks={}".format(self.total_insert_blocks)
                  )
        return output


    def dump_csv_header(self):
        output = ",".join(["total_insert_1mate_hi_qual",
                           "total_insert_1mate_lo_qual",
                           "total_insert_1mate",
                           "total_insert_agree",
                           "total_insert_agree_hi_qual",
                           "total_insert_agree_lo_qual",
                           "total_insert_conflict",
                           "total_insert_conflict_hi_qual",
                           "total_insert_conflict_hilo_qual",
                           "total_insert_conflict_lo_qual",
                           "total_insert_conflict_hino_qual",
                           "total_insert_conflict_lono_qual",
                           "total_insert_conflict_equal_hi_qual",
                           "total_inserts",
                           "total_insert_blocks",
                           "total_match_1mate",
                           "total_match_1mate_hi_qual",
                           "total_match_1mate_lo_qual",
                           "total_match_nonconflict",
                           "total_match_nonconflict_hi_qual",
                           "total_match_nonconflict_lo_qual",
                           "total_match_conflict",
                           "total_match_conflict_lo_qual",
                           "total_match_conflict_hi_qual",
                           "total_match_conflict_equal_hi_qual"])
        return output


    def dump_csv(self):
        output =  ",".join(str(x) for x in [self.total_insert_1mate_hi_qual,
                                      self.total_insert_1mate_lo_qual,
                                      self.total_insert_1mate,
                                      self.total_insert_agree,
                                      self.total_insert_agree_hi_qual,
                                      self.total_insert_agree_lo_qual,
                                      self.total_insert_conflict,
                                      self.total_insert_conflict_hi_qual,
                                      self.total_insert_conflict_hilo_qual,
                                      self.total_insert_conflict_lo_qual,
                                      self.total_insert_conflict_hino_qual,
                                      self.total_insert_conflict_lono_qual,
                                      self.total_insert_conflict_equal_hi_qual,
                                      self.total_inserts,
                                      self.total_insert_blocks,
                                      self.total_match_1mate,
                                      self.total_match_1mate_hi_qual,
                                      self.total_match_1mate_lo_qual,
                                      self.total_match_nonconflict,
                                      self.total_match_nonconflict_hi_qual,
                                      self.total_match_nonconflict_lo_qual,
                                      self.total_match_conflict,
                                      self.total_match_conflict_lo_qual,
                                      self.total_match_conflict_hi_qual,
                                      self.total_match_conflict_equal_hi_qual])
        return output