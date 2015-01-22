#! /usr/bin/env python

# Simulates an ancestral selection graph for (serial) sample data.

# February 22, 2012:
# changed because the branching was broken before.  Nodes that
# branch shouldn't have two parents with real edges to each; they
# should have one edge that forks so that the same mutations
# happen along the path to each parent.
#
# We could either have phantom nodes to represent the same lineage
# that then have two length-0 edges to the two possible parents,
# or we could have forked edges which connect the child to two parents.
# We go with the latter for now.

# February 29, 2012:
# Over the last week, changed so that annotations are now being spat
# out.  We now have to change the recursion in FigTree and
# __FT_helper__ so that it seems a little conceptually cleaner (though
# the previous version worked fine).  More to the point, we need to do
# it this way to get the labels to look right (the annotations are a
# bit of a pain, and building everything into the labels is going to
# work better in FigTree).

# March 1, 2012:
# on Art's suggestion, output the ASG as PhyloXML.  We will break cycles
# by noting that any time a not-real node is the parent of a child with
# two nodes, we break the connection; if this child has *two* non-real
# parents, then we retain whichever one ought to have been its real parent
# (based on the types).  Any broken link is noted with a CladeRelation.

# March 8, 2012:
# change the time convention: more recent samples have *bigger*
# numbers, and times are now interpreted as "days after the epoch".
# Labels must include this number, as our tree reconstruction
# software will probably need that.  Also, allow for zero selection and
# mutation rates.

# March 14, 2012:
# Add a resampling rate parameter so that the coalescence times are now
# exp(\alpha * n(n-1)/2)

# March 27, 2012:
# Changed FigTree annotations so that it includes whether or not a node
# is a leaf

# April 13, 2012:
# Tightening up the integration in calculating A1prob
# We're getting negative times, as well, so we write those as "minus"

# April 18, 2012:
# - Realized the FigTree annotations are wrong; 'taxa' should be the first column
# and it should match up with the labels in the tree.  Change what we currently
# call 'taxa' to 'nodelabel' (this is often a string representation of a tuple).
# - Modified so that it can also produce straight coalescent trees, with no
# selection


import string;
import random;
import math, scipy.integrate, scipy.stats;
import collections;
from operator import itemgetter;
#import pygraphviz as pgv;
# For PhyloXML stuff
import xml.etree.ElementTree as ET;
import numpy as np


# Nodes in the graph
class Node:
    # Each node should have:
    # - a label: tips will be labelled by numbers (or whatever), and
    #   an internal node will be labelled by a tuple containing all of the
    #   labels of its descendants
    # - allele
    # - time
    # - Edges to children
    # - Edges to parents, plural (there could be branching!)
    # Later we will introduce
    # - is_real: tells us whether this node is actually part of the tree
    #   or not

    def __init__(self, label, time):
        self.label = label;
        self.time = time;

        # This stuff is introduced in the building of the ASG
        self.child1 = None;
        self.child2 = None;
        self.parent = None;

        # This stuff is decided after the ASG is built and we trace back
        # through
        self.allele = None;

# Edges, which we use to connect nodes in coalescence events
class Edge:
    # Each edge has:
    # - the child node
    # - the parent node
    # - the length (for convenience)
    # - the number of mutations along it (decided after the ASG is built)

    def __init__(self, child, parent, num_mut = None):
        self.child = child;
        self.parent = parent;
        # March 8, 2012: self.length changed to reflect the change
        # in time convention (bigger <=> more recent)
        self.length = child.time - parent.time;
        self.num_mut = num_mut;

# Branching edges, which connect nodes in branching events
class BranchingEdge:
    # Each edge has:
    # - the child node
    # - the continuing parent
    # - the incoming parent
    # - the length (for convenience)
    # - the number of mutations along it
    # - the real parent (which is decided when both parents' types
    #   are determined)

    # PRE: c_parent and i_parent must have the same time
    def __init__(self, child, c_parent, i_parent, num_mut = None):
        self.child = child;
        self.c_parent = c_parent;
        self.i_parent = i_parent;
        # March 8, 2012: changed as Edge was
        self.length = child.time - c_parent.time;
        self.num_mut = num_mut;
        self.realparent = None;

class ASGerr:
    def __init__(self, reason):
        self.reason = reason;

class ASG:
    # This object must have:
    # - all of the currently active nodes
    # - a current time
    # - a list of when new samples should be introduced, and how many
    #   * this can be stored as a deque of tuples (when, how many)
    #     and should be sorted by when in ascending order so that
    #     the first entry is the first thing to be introduced and so on.
    #     This is more efficient than using a list as we want to pop
    #     from the head of the deque and append to the tail.
    # - a list (stack) of all nodes in the tree, so that we can
    #   traverse the tree from top to bottom

    # Helper functions: these are used to compute the density of the
    # stationary distribution of X_1(t), which is the diffusion
    # approximation of the fraction of A_1 types in the population
    # and is used in determining the type of the ultimate ancestor
    def __unnormalized_density__(self, x):
        poly_factor = pow(x*(1-x), self.mutation_rate - 1);
        return poly_factor * math.exp(-self.selection_rate * x);
    
    def stat_density(self, x):
        return self.__unnormalized_density__(x)/self.normalizing_const;
    

    # whens and how_manys are lists that should match up
    def __init__(self, resample_rate, selection_rate, mutation_rate,
                 whens, how_manys, growth_rate = None):
        self.curr_lineages = collections.OrderedDict()  # use OrderedDict to ensure reproducability when random_seed set
        # March 8, 2012: set time to be the most recent time in whens;
        # i.e. the largest one
        self.time = max(whens);
        self.last_label = 0;
        self.nodes_sorted = [];

        # There must be at least one sample
        if len(whens) < 1:
            raise ASGerr("There must be at least one sample specified");

        # Make a list of tuples (when, how_many)
        if len(whens) != len(how_manys):
            raise ASGerr("Lists of times and sizes of samples do not have the same length");

        whens_how_manys = [];
        for i in range(len(whens)):
            # Check that all samples are of size > 0;
            if how_manys[i] <= 0:
                raise ASGerr("Samples have non-positive sizes");
            whens_how_manys += [(whens[i], how_manys[i])];
        # March 8, 2012: now we take larger times to be more recent
        whens_how_manys = sorted(whens_how_manys, key=itemgetter(0),
                                 reverse = True);
        
        self.introduce_new = collections.deque(whens_how_manys);
        if selection_rate < 0:
            raise ASGerr("Negative selection rate");
        self.selection_rate = selection_rate;
        if mutation_rate < 0:
            raise ASGerr("Negative mutation rate");
        self.mutation_rate = mutation_rate;
        if resample_rate < 0:
            raise ASGerr("Negative resampling rate");
        self.resample_rate = resample_rate;
        
        
        if growth_rate is not None and growth_rate < 0:
            raise ASGerr("Negative population growth rate");
        
        # logistic growth rate parameters
        self.growth_rate = growth_rate
        self.carrying_cap = 1./resample_rate

        # Figure out what the normalizing constant must be.  We do
        # this here, exactly once, as this doesn't need to be
        # recomputed every time you call stat_density.
        # -- RL 2012_04_13
        self.normalizing_const = 1;
        self.A1prob = 0.5;
        if selection_rate != 0:
            self.normalizing_const = scipy.integrate.quad(self.__unnormalized_density__, 0, 1)[0];

            # Find the probability of the ultimate ancestor being of type A1.
            self.A1prob = scipy.integrate.quad(lambda y: y*self.stat_density(y),
                                               0, 1)[0];
        print("Probability of UA being A1: " + str(self.A1prob));


    # Builds the ASG to the next event.  Returns True if we have reached
    # the ultimate ancestor.
    def next_step(self):

        # Generate random waiting times for the next branching and coalescence
        # event
        num_lineages = len(self.curr_lineages);

        # Quick check in case we call this after the ASG is complete.
        # In this case we do nothing.
        if num_lineages == 1 and len(self.introduce_new) == 0:
            return True;

        # If this is zero, then we know we have to introduce the first
        # tip samples to the system.  By construction the first samples
        # must be added at time 0.
        if num_lineages == 0:
            when, how_many = self.introduce_new.popleft();
            
            # Add these first nodes to curr_lineages, with no children
            # and no parents.
            for i in range(how_many):
                # Label it with self.last_label + 1
                self.curr_lineages[self.last_label + 1] = Node(self.last_label+1, when);
                self.nodes_sorted.append(self.curr_lineages[self.last_label+1]);
                self.last_label += 1;
                
            # Trivial case: if there is only one sample of size 1, then we
            # are done
            if len(self.curr_lineages) == 1 and len(self.introduce_new) == 0:
                return True;
            return False;

        branch_wait = -1;
        if self.selection_rate > 0:
            branch_wait = random.expovariate(float(num_lineages)*self.selection_rate/2.0);

        coal_wait = -1;
        # If there's only one lineage left but more to come, there can
        # be no coalescence.  We leave coal_wait set at -1 in that case.
        if num_lineages > 1:
            if self.growth_rate is None:
                coal_wait = random.expovariate(self.resample_rate * float(num_lineages*(num_lineages-1))/2.0);
            else:
                # implicitly assumes population size starts at 1
                pop_size = self.carrying_cap * math.exp(self.growth_rate * self.time)
                pop_size /= self.carrying_cap + math.exp(self.growth_rate * self.time) - 1
                if num_lineages % 100 == 0:
                    print self.time, pop_size
                
                coal_wait = random.expovariate(1./pop_size * float(num_lineages*(num_lineages-1))/2.0);

        # Case 1: there is branching and more than one lineage.
        if self.selection_rate > 0 and num_lineages > 1:
            curr_wait = min(branch_wait, coal_wait);
        # Case 2: there is branching and only one lineage
        elif self.selection_rate > 0 and num_lineages == 1:
            curr_wait = branch_wait;
            # Set coal_wait to be bigger than branch_wait
            coal_wait = branch_wait + 1;
        # Case 3: no branching, more than one lineage.
        elif self.selection_rate == 0 and num_lineages > 1:
            curr_wait = coal_wait;
            # Set branch_wait to be bigger than coal_wait
            branch_wait = coal_wait + 1;
        # Case 4: no branching, only one lineage
        else:
            # then we set curr_wait to be big enough so that it
            # introduces new samples
            curr_wait = self.time - self.introduce_new[0][0] + 1;

        # Check if this will make the time exceed the next time of
        # samples, if there are more samples
        if len(self.introduce_new) > 0:
            if self.time - curr_wait < self.introduce_new[0][0]:
                # If so, then rather than do any coalescence or branching,
                # we introduce new nodes to the system, essentially just
                # like we do in introducing the first samples.
                when, how_many = self.introduce_new.popleft();

                for i in range(how_many):
                    # Label it with self.last_label + 1
                    self.curr_lineages[self.last_label + 1] = Node(self.last_label+1, when);
                    self.nodes_sorted.append(self.curr_lineages[self.last_label+1]);
                    self.last_label += 1;
                # Increment the time
                self.time = when;

                # This should never happen, but check to see if we are out of
                # tips to add, and that we are down to one lineage
                if len(self.curr_lineages) == 1 and len(self.introduce_new) == 0:
                    raise ASGerr("Only one lineage after adding more tips");
            
                return False;

        # Decide whether this is a coalescence or a branching.
        return_val = False;
        if coal_wait < branch_wait:
            # Choose two lineages to coalesce
            coal_labels = random.sample(self.curr_lineages.keys(), 2);

            # Combine their labels
            new_label_list = [];
            for label in coal_labels:
                if type(label) == tuple:
                    for entry in label:
                        new_label_list += [entry];
                else:
                    new_label_list += [label];
            new_label = tuple(new_label_list);

            new_node = Node(new_label, self.time - curr_wait);
            self.curr_lineages[new_label] = new_node;
            self.nodes_sorted.append(new_node);

            child1_edge = Edge(self.curr_lineages[coal_labels[0]],
                               new_node);
            child2_edge = Edge(self.curr_lineages[coal_labels[1]],
                               new_node);
            
            new_node.child1 = child1_edge;
            new_node.child2 = child2_edge;

            self.curr_lineages[coal_labels[0]].parent = child1_edge;
            self.curr_lineages[coal_labels[1]].parent = child2_edge;

            # Now remove the old coalesced nodes from curr_lineages
            # (they'll still be referred to by their parent);
            self.curr_lineages.pop(coal_labels[0]);
            self.curr_lineages.pop(coal_labels[1]);

            # Check that we've reached the UA
            if len(self.curr_lineages) == 1 and len(self.introduce_new) == 0:
                return_val = True;

        else:
            # Choose a lineage to branch
            coal_label = random.choice(self.curr_lineages.keys());
            branching_node = self.curr_lineages[coal_label];

            # Make the parents
            c_parent = Node(self.last_label+1, self.time - curr_wait);
            i_parent = Node(self.last_label+2, self.time - curr_wait);
            self.last_label += 2;

            branching_edge = BranchingEdge(branching_node, c_parent,
                                           i_parent);

            # Introduce children and parents to each other
            branching_node.parent = branching_edge;
            c_parent.child1 = branching_edge;
            i_parent.child1 = branching_edge;

            # Add parents to curr_lineages, remove child
            self.curr_lineages[c_parent.label] = c_parent;
            self.curr_lineages[i_parent.label] = i_parent;
            self.curr_lineages.pop(coal_label);
            self.nodes_sorted.append(c_parent);
            self.nodes_sorted.append(i_parent);

            # Sanity check: this should never happen
            if len(self.curr_lineages) == 1 and len(self.introduce_new) == 0:
                raise ASGerr("Reached UA on a branching event");

        # Decrement the time
        self.time -= curr_wait;
        return return_val;

    # Function to build the ASG up to the time it reaches the ultimate
    # ancestor.
    #
    # If random_seed is specified then the simulation will use this seed
    # (this is mostly for debugging purposes).
    def buildASG(self, random_seed=None):
        random.seed(random_seed);
        np.random.seed(random_seed)  # Need this to reproduce random variables generated by scipy

        while not self.next_step():
            pass

    
    # This function will assign types to the tree.  There are two types,
    # types A1 and A2, with A2 being advantageous.
    # PRE: the ASG must be fully built
    def assign_types(self):
        decide_root = random.uniform(0,1);

        root_type = "A1";
        if decide_root > self.A1prob:
            root_type = "A2";

        # Get the root
        root_node = self.curr_lineages[self.curr_lineages.keys()[0]];
        root_node.allele = root_type;

        # Now go through the tree, top-to-bottom, assigning types.
        # This basically means going through nodes_sorted backwards,
        # skipping the last entry (that's the root).
        for i in range(len(self.nodes_sorted)-2, -1, -1):
            curr_node = self.nodes_sorted[i];
            p_edge = curr_node.parent;

            p_edge.num_mut = 0;
            if self.mutation_rate > 0:
                p_edge.num_mut = scipy.stats.poisson.rvs(p_edge.length * self.mutation_rate);
            
            # Case 1: this node has one parent (so it's part of
            # a coalescence event going back in time)
            if p_edge.__class__.__name__ == "Edge":
                parent = p_edge.parent;
                # Decide the type
                if p_edge.num_mut % 2 == 0:
                    curr_node.allele = parent.allele;
                else:
                    # Flip the type
                    if parent.allele == "A1":
                        curr_node.allele = "A2";
                    else:
                        curr_node.allele = "A1";
                        
            # Case 2: this node has two parents (so it's part of a
            # branching event going back in time)
            else:
                c_parent = p_edge.c_parent;
                i_parent = p_edge.i_parent;

                c_type = c_parent.allele;
                i_type = i_parent.allele;

                # Now: with this information we can decide for the
                # current node what its allele should be and which
                # parent is its real parent.
                #
                # ancestral_type is the type of the lineage right
                # after (i.e. moving forwards in time, so right below
                # the branch point) branching.
                ancestral_type = "A1";
                if c_type == "A2" or i_type == "A2":
                    ancestral_type = "A2";

                if p_edge.num_mut%2 == 0:
                    curr_node.allele = ancestral_type;
                else:
                    if ancestral_type == "A1":
                        curr_node.allele = "A2";
                    else:
                        curr_node.allele = "A1";

                p_edge.realparent = c_parent;
                if i_type == "A2":
                    p_edge.realparent = i_parent;


    # This routine goes through and marks the true genealogy of the
    # samples.
    #
    # A node is included if:
    #  - it is a leaf
    #  - it is a coalescence node and at least one of its children is real
    #  - it is a branching node, is the real parent of its child, and its
    #    child is real
    # 
    # PRE: ASG must be built and types must be assigned
    def mark_true_genealogy(self):
        # Go through the list of nodes, from bottom to top, and mark
        # things as real or not real.

        for node in self.nodes_sorted:
            if node.child1 == None and node.child2 == None:
                # This is a leaf
                node.is_real = True;
            elif node.child1.__class__.__name__ == "Edge":
                # This is a coalescence node; it's real if either of
                # its children are real
                node.is_real = node.child1.child.is_real or node.child2.child.is_real;
            elif node.child1.__class__.__name__ == "BranchingEdge":
                # This is a branching node; it's real if it is the real
                # parent of its child
                node.is_real = node.child1.realparent == node and node.child1.child.is_real;


    # Now, a routine that outputs the true genealogy in Newick format
    # (with FigTree extensions).
    def FigTree(self):
        # We will need a recursive helper.
        root_node = self.nodes_sorted[-1];
        return self.__FT_helper__(root_node, 0, 0);

    # The recursive helper.  Algorithm:
    # - any time you go left, output a left parens
    # - any time you return back up a branch, output a right parens,
    #   and then a colon, comments, and a length
    # - when you reach a leaf, output the letter
    #
    # dist_to_parent, mut_to_parent are accumulator variables that
    # track how far base_node is from its last "relevant" ancestor in
    # the tree.
    #
    # A node is bifurcating if it has two children, both real.
    # This returns a 2-tuple: the first element is the Newick string
    # representing the tree rooted at this node (minus the length);
    # the second is a list of annotations for every node below this one,
    # inclusive
    def __FT_helper__(self, base_node, dist_to_parent, mut_to_parent):

        # DEBUGGING:
        #print("Calling FT_helper on node " + self.__Newick_label__(base_node.label));

        # Base case: if base_node is a leaf, output the label and
        # the annotation for this node
        if base_node.child1 == None and base_node.child2 == None:
            leaf_label = self.__Newick_label__(base_node.label,
                                               base_node.allele,
                                               mut_to_parent,
                                               base_node.time);
            leaf_annot = [(leaf_label, True, base_node.label, base_node.allele,
                            mut_to_parent, base_node.time)];
            leaf_Newick_temp = string.Template("$label:$dist");
            return (leaf_annot,
                    leaf_Newick_temp.substitute(label=leaf_label,
                                                 dist=dist_to_parent));

        # Base case 2: if this node is not real, raise an exception.
        # This should never actually happen.
        if not base_node.is_real:
            raise ASGerr("__FT_helper__ was called on a non-real node");

        # Now: we know that this is a real internal node, so we
        # recurse down to the children.  We can't recurse down right
        # away before testing whether the children are real; we need
        # to know if this current node is a relevant ancestor to its
        # children or not.
        # The following are possible:
        # - there are two real children; in which case, this is a
        #   branch point and we recurse down to both
        # - there is only one real child, in which case we recurse down
        #   on that child
        firstchild = base_node.child1.child;

        sc_is_real = False;
        secondchild = None;
        if base_node.child2 != None:
            secondchild = base_node.child2.child;
            sc_is_real = secondchild.is_real;

        if firstchild.is_real and sc_is_real:
            # Recurse down to both, resetting the values of dist_to_parent
            # in the recursive call
            
            fc_annot, fc_tree = self.__FT_helper__(firstchild,
                                                   firstchild.parent.length,
                                                   firstchild.parent.num_mut);
            sc_annot, sc_tree = self.__FT_helper__(secondchild,
                                                   secondchild.parent.length,
                                                   secondchild.parent.num_mut);

            base_Newick_label = self.__Newick_label__(base_node.label,
                                                      base_node.allele,
                                                      mut_to_parent,
                                                      base_node.time);
            base_annot = (base_Newick_label, False, base_node.label,
                          base_node.allele, mut_to_parent, base_node.time);
            all_annot = fc_annot + sc_annot + [base_annot];
            subtree_template = string.Template("($firstchild,$secondchild)$base_label:$base_length");

            return (all_annot,
                    subtree_template.substitute(firstchild=fc_tree,
                                                 secondchild=sc_tree,
                                                 base_label=base_Newick_label,
                                                 base_length=dist_to_parent));

        # Last case: only one child is real.  Recurse on that child,
        # incrementing the values of dist_to_parent and mut_to_parent
        realchild = firstchild;
        if sc_is_real:
            realchild = secondchild;

        added_dist = realchild.parent.length;
        added_mut = realchild.parent.num_mut;
        return self.__FT_helper__(realchild,
                                  dist_to_parent + added_dist,
                                  mut_to_parent + added_mut);


    # Helper that translates a label (which may be a tuple) into a
    # valid Newick name; it turns (5,2), allele, mut, time into
    # 5_2__type[allele]_mut[mut]_time[time]
    # Replace minus signs on negative times with "neg" -- 2012_04_13
    def __Newick_label__(self, label, allele, mut_to_parent, time):
        newick_label = "";
        if type(label) == tuple:
            newick_label +=  str(label[0]);
            for entry in label[1:]:
                newick_label += "_" + str(entry);
        else:
            newick_label = str(label);

        newick_label += "__type" + allele;
        newick_label += "_mut" + str(mut_to_parent);
        newick_label += "_time" + str(time).replace(".", "pt").replace("-", "neg");

        return newick_label;

    # Another recursive function that will output a representation of the
    # entire ASG as PhyloXML.  Cycles are broken by picking the "real"
    # of the two parents any time a child has two parents; the relationship
    # between the "fake" parent and the child is noted as a CladeRelation.
    def PhyloXML(self):

        phylo_root = ET.Element("Phylogeny");

        # Call a recursive helper on the UA of the ASG to do the rest.
        # The clade relations will come as a list, which we will put
        # together here
        ult_anc = self.nodes_sorted[-1];
        ancestral_tree, clade_relations = self.__PXML_helper__(ult_anc);
        
        phylo_root.append(ancestral_tree);
        phylo_root.extend(clade_relations);

        phylo_tree = ET.ElementTree(phylo_root);
        return phylo_tree;

    # The workhorse function for the above.  Returns a tuple with
    # the tree rooted at this node, as well as a list of CladeRelations
    # describing the "fake parent-child" relationships under this mode
    def __PXML_helper__(self, base_node):

        # Describe this current clade
        curr_clade = ET.Element("Clade");

        # If this is the root, set branch_length = 0; otherwise
        # set it to the branch length
        if base_node.parent == None:
            curr_clade.set("branch_length", "0");
        else:
            curr_clade.set("branch_length", str(base_node.parent.length));
        curr_clade.set("id_source", str(base_node.label));
        curr_name = ET.Element("Name");
        curr_name.text = str(base_node.label);
        curr_clade.append(curr_name);

        # If this is not the root, add a property for the number of
        # mutations on the parent branch.
        if base_node.parent != None:
            curr_mut = ET.Element("Property");
            curr_mut.set("ref", "mutations");
            curr_mut.set("datatype", "xsd:integer");
            curr_mut.set("applies_to", "parent_branch");
            curr_mut.text = str(base_node.parent.num_mut);
            curr_clade.append(curr_mut);
        
        curr_allele = ET.Element("Property");
        curr_allele.set("ref", "allele");
        curr_allele.set("datatype", "xsd:string");
        curr_allele.set("applies_to", "node");
        curr_allele.text = base_node.allele;
        
        curr_is_real = ET.Element("Property");
        curr_is_real.set("ref", "is_real");
        curr_is_real.set("datatype", "xsd:boolean");
        curr_is_real.set("applies_to", "node");
        curr_is_real.text = str(base_node.is_real);
        
        curr_time = ET.Element("Date");
        curr_time.set("unit", "days after epoch");
        time_value = ET.SubElement(curr_time, "value");
        time_value.text = str(base_node.time);
        
        curr_clade.extend([curr_allele, curr_is_real, curr_time]);

        # Base case: this is a leaf.  In this case, return a single clade
        # representing this node, and an empty list of clade relations.
        # FIXME: we may need to stringify some fields later
        if base_node.child1 == None and base_node.child2 == None:
            return (curr_clade, []);

        # Recursive case: this is an internal node.  This is either a
        # coalescence node or the result of a branching node; we check
        # by the number of children.
        
        curr_clade_rels = [];
        
        firstchild = base_node.child1.child;
        secondchild = None;
        if base_node.child2 != None:
            secondchild = base_node.child2.child;

        if secondchild == None:
            # This is the result of a branching node.  We check to see if
            # it is the real parent of its child.
            if firstchild.parent.realparent == base_node:
                # Recurse on the child and add it to curr_clade.
                # There is no false parent relationship here.
                child_subtree, child_clade_rels = self.__PXML_helper__(firstchild);
                curr_clade.append(child_subtree);
                curr_clade_rels = child_clade_rels;                

            else:
                # Return curr_clade as well as a CladeRelation
                # indicating the false parent-child relationship
                f_p_r_attrib = {"id_ref_0": str(base_node.label),
                                "id_ref_1": str(firstchild.label),
                                "type": "false parent"};
                false_parent_r = ET.Element("CladeRelation",
                                            attrib = f_p_r_attrib);
                curr_clade_rels = [false_parent_r];
                
        else:
            # This is a coalescence node.  Recurse on the children and
            # append the false parent relationships from both.
            child1_subtree, child1_clade_rels = self.__PXML_helper__(firstchild);
            child2_subtree, child2_clade_rels = self.__PXML_helper__(secondchild);
            curr_clade.extend([child1_subtree, child2_subtree]);
            curr_clade_rels = child1_clade_rels + child2_clade_rels;

        return (curr_clade, curr_clade_rels);

    # Function to make NeXML output, which is more appropriate for the
    # ASG.  It's more accommodating to graph structure and is laid out
    # more like Graphviz input.
    def NeXML(self):
        NeXMLroot = ET.Element("nex:nexml");
        NeXMLroot.set("version", "0.9");
        NeXMLroot.set("generator", "ASGsim.py");
        NeXMLroot.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
        NeXMLroot.set("xmlns:xml", "http://www.w3.org/XML/1998/namespace");
        NeXMLroot.set("xmlns:nex", "http://www.nexml.org/2009");
        NeXMLroot.set("xmlns", "http://www.nexml.org/2009");
        NeXMLroot.set("xsi:schemaLocation", "http://www.nexml.org/2009 http://www.nexml.org/2009/nexml.xsd");

        edges = [];
        nodes = [];

        # Go through nodes_sorted and add edges and nodes from there
        for curr_node in self.nodes_sorted:
            new_XMLnode = ET.Element("node");
            new_XMLnode.set("id", str(curr_node.label));
            new_XMLnode.set("label", str(curr_node.label));
            if curr_node.parent == None:
                # This is the root
                new_XMLnode.set("root", "true");

            # Add meta subelements for allele, is_real
            allele_meta = ET.Element("meta");
            allele_meta.set("property", "allele");
            allele_meta.set("datatype", "xsd:string");
            allele_meta.set("xsi:type", "nex:LiteralMeta");
            allele_meta.set("content", curr_node.allele);
            new_XMLnode.append(allele_meta);

            is_real_attrib = {"property": "is_real",
                              "datatype": "xsd:boolean",
                              "content": str(curr_node.is_real),
                              "xsi:type": "nex:LiteralMeta"};
            ET.SubElement(new_XMLnode, "meta", attrib=is_real_attrib);

            time_attrib = {"property": "time",
                           "datatype": "xsd:double",
                           "content": str(curr_node.time),
                           "xsi:type": "nex:LiteralMeta"};
            ET.SubElement(new_XMLnode, "meta", attrib=time_attrib);
            
            nodes += [new_XMLnode];

            # If this node has children, add the edges to the list of
            # edges.  "source" is the parent and "target" is the child.
            for curr_edge in [curr_node.child1, curr_node.child2]:
                if curr_edge != None:
                    new_XMLedge = ET.Element("edge");
                    new_XMLedge.set("id", str(curr_node.label) + "_" + str(curr_edge.child.label));
                    new_XMLedge.set("source", str(curr_node.label));
                    new_XMLedge.set("target", str(curr_edge.child.label));
                    new_XMLedge.set("length", str(curr_edge.length));

                    # Add meta subelements for number of mutations,
                    # and whether this edge connects the real parent to the
                    # child or the fake parent
                    num_mut_attrib = {"property": "mutations",
                                      "content": str(curr_edge.num_mut),
                                      "datatype": "xsd:integer",
                                      "xsi:type": "nex:LiteralMeta"};
                    ET.SubElement(new_XMLedge, "meta", attrib=num_mut_attrib);

                    is_realparent = True;
                    if curr_edge.__class__.__name__ == "BranchingEdge":
                        is_realparent = curr_edge.realparent == curr_node;
                    realparent_attr = {"property": "realparent",
                                       "content": str(is_realparent),
                                       "datatype": "xsd:boolean",
                                       "xsi:type": "nex:LiteralMeta"};
                    ET.SubElement(new_XMLedge, "meta", attrib=realparent_attr);

                    edges += [new_XMLedge];
            
            
        # Now make a network and add the nodes and edges.
        XML_ASG = ET.Element("network",
                             attrib = {"id": "ASG",
                                        "xsi:type": "nex:FloatNetwork"});
        XML_ASG.extend(nodes + edges);
        XML_trees = ET.Element("trees", attrib = {"id": "ASGs"});
        XML_trees.append(XML_ASG);
        NeXMLroot.append(XML_trees);

        return ET.ElementTree(NeXMLroot);
