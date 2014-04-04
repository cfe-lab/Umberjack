import slice_miseq
import sam_handler
import os
import subprocess
import HyPhy


def process_windows(ref_fasta_filename, sam_filename, mapping_cutoff, read_qual_cutoff, max_prop_N, cov_thresh):
    """
    Slides window along genome.
    Creates the multiple sequence aligned fasta files for the window.
    Feeds the multiple-sequence aligned fasta files to fasttree2 to create a tree.
    Feeds the tree into HyPhy to obtain dn/ds values.
    """
    # TODO:  how to get ORFs?
    # TODO:  do not use kalign2 from PATH
    # TODO:  do not hyphy from PATH

    # keep track of reference contig/chromosome and its length
    ref_to_len = {}
    with open(ref_fasta_filename, 'r') as ref_fasta_fh:
        for line in ref_fasta_fh:
            line = line.rstrip().split()[0]
            header = ''
            seq_len = 0
            if line[0] == '>':
                if seq_len:
                    ref_to_len[header] = seq_len

                header = line[1:]
            else:
                seq_len += len(line)

        ref_to_len[header] = seq_len

    # For every contig/chromosome in the reference fasta, create a pseudo multiple-sequence aligned fasta file
    # TODO:  handle indels in multiple sequence align file
    # TODO:  allow option to specify msa file location
    sam_filename_prefix, fileExtension = os.path.splitext(sam_filename)
    msa_fasta_filename = sam_filename_prefix + ".msa.fasta"

    # TODO:  need to create separate msa fasta file for each reference contig/chromosome
    sam_handler.get_msa_fasta_from_sam(sam_filename=sam_filename, mapping_cutoff=mapping_cutoff,
                                           read_qual_cutoff=read_qual_cutoff, max_prop_N=max_prop_N,
                                           out_fasta_filename=msa_fasta_filename)

    # Slice the multiple sequence aligned fasta file and feed to fasttree2
    # TODO:  handle multiple contig/chromosome in reference fasta.  Right now only handles one contig/chromosome.
    windowsize = slice_miseq.get_best_window_size_by_read_cov(sam_filename=sam_filename, ref_filename=ref_fasta_filename, cov_thresh=cov_thresh)

    for ref in ref_to_len:
        last_window_start_pos = ref_to_len[ref] - windowsize
        for start_pos in range(1, last_window_start_pos):  # start_pos is 1-based
            end_pos = start_pos + windowsize                # end_pos is 1-based

            fastree_logfilename = sam_filename_prefix + ".fasttree.log"
            fastree_treefilename = sam_filename_prefix + ".tree"
            fasttree_stdouterr_filename = sam_filename_prefix + ".stdouterr.txt"

            with open(fasttree_stdouterr_filename, 'w') as fasttree_stdouterr_fh:
                slice_miseq.slice_msa_fasta(fasta_filename=msa_fasta_filename, start_pos=start_pos, end_pos=end_pos)
                subprocess.check_call(['FastTree', '-gtr', '-nt', '-nosupport',
                                       '-log', fastree_logfilename, '-out', fastree_treefilename,
                                       msa_fasta_filename],
                                      stdout=fasttree_stdouterr_fh, stderr=fasttree_stdouterr_fh, shell=False)



            # pass in fasttree2 tree to HyPHy
            # SetDialogPrompt ("Select file containing paths to FASTA files");
            # fscanf (PROMPT_FOR_FILE, "Lines", paths);
            #
            # for (f = 0; f < Columns(paths); f = f + 1)
            # {
            #         fprintf (stdout, "processing ", paths[f], "\n");
            #
            #         stdinRedirect = {};
            #
            #         stdinRedirect["00"] = "Universal";
            #         stdinRedirect["01"] = paths[f];
            #         stdinRedirect["02"] = "Select from list";
            #         stdinRedirect["03"] = "V3, clinical";
            #         stdinRedirect["04"] = "HIV 25%";
            #         stdinRedirect["05"] = "All";
            #         stdinRedirect["06"] = "CSV";         /* insert reference sequence? */
            #         stdinRedirect["07"] = paths[f]^{{"\\.prealign\\.fas"}{".csv"}};
            #
            #         fprintf (stdout, stdinRedirect, "\n");
            #
            #         ExecuteAFile ("/Users/apoon/wip/emeline/emelineV2.6_SA.bf", stdinRedirect);
            # }

