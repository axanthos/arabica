from __future__ import division

"""arabica.py -- unsupervised learning of root-and-pattern morphology in Python
(c) 2016 Aris Xanthos, distributed under GNU GPL v3 license.
"""

from pytrie import StringTrie

import math
import copy
import os

__version__ = 2.0       # 1.0 was in Perl

# Parameters...
input_file         = os.path.join('../data/raw_corpus.txt')
max_num_iterations = 10
del_duplicates     = False
ins_symbol         = '-'

def main():

    # Read word list.
    word_list = read_word_list(input_file, del_duplicates)

    # Initial decomposition of roots and patterns.
    initial_RP_decomp = Initial_RP_decomp(word_list, phonological_strategy)

    # Initialize RP analysis.
    my_RP_analysis = RP_analysis(initial_RP_decomp)

    # Apply bootstrap heuristics.
    my_RP_analysis.bootstrap_heuristics(
        min_root_length=2,
        min_vowels_in_pattern=1,
        min_num_roots=2,
        min_num_patterns=2,
    )

    # Incremental heuristics...
    iter = 1
    while (iter <= max_num_iterations):

        initial_description_length = my_RP_analysis.description_length()

        # Extend known RP structures.
        my_RP_analysis.extend_known_RP_structures(iter)

        # Extend fragments of known RP structures.
        my_RP_analysis.extend_fragments_known_RP_structures(
            min_num_patterns=2,
            iter=iter,
        )

        # Extend known roots to known patterns.
        my_RP_analysis.extend_known_roots_to_known_patterns(iter)
        
        # Extend known roots.
        my_RP_analysis.extend_known_roots(
            min_robustness=5,
            min_pattern_count=2,
            iter=iter,
        )
        
        # Stop iterating if description length didn't decrease...
        if initial_description_length <= my_RP_analysis.description_length():
            break

        iter += 1

    print my_RP_analysis.display()
    


def read_word_list(input_file, del_duplicates):
    """Attempt to read word list from file and return them in list form"""
    word_list = [line.strip() for line in open(input_file)]
    if del_duplicates:
        word_list = list(set(word_list))
    return word_list


class Initial_RP_decomp(object):
    """
    Abstraction for the operation of providing an initial RP decomposition
    of a word list.
    """
    def __init__(self, word_list, RP_decomp_strategy):
        self.word_list = word_list
        self.consonants, self.vowels = learn_CV(word_list)
        self.word_decomp = RP_decomp_strategy(self)


def frequency_based_strategy():
    """RP decomposition method from Bilal & Johnson 2013"""
    pass


def phonological_strategy(initial_RP_decomp):
    """RP decomposition method based on Sukhotin's algorithm"""
    word_decomp = list()
    for word in initial_RP_decomp.word_list:
        root = ''.join([s for s in word if s in initial_RP_decomp.consonants])
        pattern = ''.join(
            [s if s in initial_RP_decomp.vowels else ins_symbol for s in word]
        )
        word_decomp.append((root, pattern))
    return word_decomp


def learn_CV(word_list):
    """Apply Sukhotin's algorithm to identify consonants and vowels"""

    print "Learning CV..."

    # Compute cooccurrences...
    cooccurrences = dict()
    symbols       = set()
    for word in word_list:
        symbols.update(word)
        for pos in xrange(len(word)-1):
            left_symbol   = word[pos]
            right_symbol  = word[pos+1]
            if left_symbol != right_symbol:
                try:
                    cooccurrences[(left_symbol, right_symbol)] += 1
                    cooccurrences[(right_symbol, left_symbol)] += 1
                except KeyError:
                    cooccurrences[(left_symbol, right_symbol)]  = 1
                    cooccurrences[(right_symbol, left_symbol)]  = 1
    symbols = list(symbols)

    # Initialize the is_consonant dict to True.
    is_consonant = dict((s, True) for s in symbols)

    # Initialize the vocalic_score dict to total cooccurrence frequency...
    vocalic_score = dict((s, 0) for s in symbols)
    for symbol in symbols:
        for symbol2 in symbols:
            try:
                vocalic_score[symbol] += cooccurrences[(symbol, symbol2)]
            except KeyError:
                pass

    # Get consonant with highest vocalic score...
    candidate = max(
        [k for (k, v) in is_consonant.iteritems() if v],
        key=vocalic_score.get
    )
    max_score = vocalic_score[candidate]

    # Iterative detection of vowels...
    while max_score > 0:

        # Re-classify candidate as vowel.
        is_consonant[candidate] = False

        # Re-compute vocalic scores of remaining consonants accordingly...
        for symbol in symbols:
            if is_consonant[symbol]:
                try:
                    vocalic_score[symbol] -= (
                        2 * cooccurrences[(symbol, candidate)]
                    )
                except KeyError:
                    pass

        # Get consonant with highest vocalic score...
        candidate = max(
            [s for (s, v) in is_consonant.iteritems() if v],
            key=vocalic_score.get
        )
        max_score = vocalic_score[candidate]

    # Build and return lists of consonants and vowels...
    consonants = list()
    vowels     = list()
    for symbol in sorted(symbols, key=lambda s: vocalic_score[s], reverse=True):
        if is_consonant[symbol]:
            consonants.append(symbol)
        else:
            vowels.append(symbol)
    return consonants, vowels

def check_RP_decomp(
        RP_decomp,
        min_stem_length=1,
        min_num_root_stems=1,
        min_num_root_suffixes=1,
):
    """Try to transfer material between roots and patterns (using successor
    count).
    
    This function is work in progress and is not currently run by default.
    """
    print "Checking RP decomposition..."

    # Initializations...
    bases = list()
    suffixes = list()
    root_stem = dict()
    root_suffix = dict()
    root_analysis = list()
    suffixes_for_stem = dict()
    stems_for_suffixes = dict()

    # Get list of consonantal skeletons and store it in trie...
    roots = [r for (r, p) in RP_decomp.word_decomp]
    trie = StringTrie(
        zip(roots, [0 for root in roots])
    )

    # For each root token...
    for root in roots:

        # If root has known segmentation, use it...
        if root in root_stem:
            root_analysis.append((root_stem[root], root_suffix[root]))
            continue

        # Else apply successor count algorithm...
        successor_counts = list()
        for length in xrange(1, len(root)+1):
            substring = root[:length]
            successor_counts.append(
                len(
                    set(
                        s[length] if len(s) > len(substring) else '#'
                        for s in trie.keys(substring)
                    )
                )
            )
        pos = len(root)-2
        while (pos > min_stem_length - 2):
            if (
                successor_counts[pos] > 1 and
                successor_counts[pos] >= successor_counts[pos-1] and
                successor_counts[pos-1] < 3 and
                successor_counts[pos+1] == 1
            ):
                root_stem[root] = root[0:pos+1]
                root_suffix[root] = root[pos+1:]
                root_analysis.append((root[0:pos+1], root[pos+1:]))
                break
            pos -= 1

        # Default is no segmentation...
        if root not in root_stem:
            root_stem[root] = root
            root_suffix[root] = ''
            root_analysis.append((root, ''))

    # Build signatures...
    for stem, suffix in root_analysis:
        try:
            suffixes_for_stem[stem].append(suffix)
        except KeyError:
            suffixes_for_stem[stem] = [suffix]
    for stem in suffixes_for_stem:
        suffixes = tuple(sorted(set(suffixes_for_stem[stem])))
        try:
            stems_for_suffixes[suffixes].append(stem)
        except KeyError:
            stems_for_suffixes[suffixes] = [stem]

    # Remove signatures with insufficiently many stems and suffixes...
    suffixes_to_delete = list()
    for suffixes in stems_for_suffixes:
        if (
            len(suffixes) < min_num_root_suffixes or
            len(stems_for_suffixes[suffixes]) < min_num_root_stems
        ):
            suffixes_to_delete.append(suffixes)
    for suffixes in suffixes_to_delete:
        del stems_for_suffixes[suffixes]

    # Build replacement mapping...
    replace_root = dict()
    for suffixes in stems_for_suffixes:
        for suffix in suffixes:
            if suffix == '':
                continue
            for stem in stems_for_suffixes[suffixes]:
                replace_root[stem + suffix] = (stem, suffix)

    # Transferring suffix to pattern for each word...
    for word_idx in xrange(len(RP_decomp.word_list)):
        word_root, word_pattern = RP_decomp.word_decomp[word_idx]
        if word_root in replace_root:
            stem, suffix = replace_root[word_root]
            replacement_root = stem
            replacement_pattern = word_pattern
            for i in xrange(len(suffix)):
                last_insertion_point = replacement_pattern.rfind(ins_symbol)
                tmp_list = list(replacement_pattern)
                tmp_list[last_insertion_point] = suffix[-(i+1)]
                replacement_pattern = ''.join(tmp_list)
            RP_decomp.word_decomp[word_idx] = (
                replacement_root,
                replacement_pattern,
            )


class RP_analysis(object):
    """Abstraction for a (model, data description) pair"""
    def __init__(self, initial_RP_decomp, initialize=True):
        self.initial_RP_decomp = initial_RP_decomp
        self.roots = list()
        self.patterns = list()
        self.RP_structures = list()
        self.data_description = list()
        self.null_model_length = None
        if initialize:
            self.initialize()

    def initialize(self):
        """Initialize RP analysis with null model"""
        creator_tag = 'INIT'
        print "Initializing RP analysis..."

        # Create null pattern...
        null_pattern = Pattern(self.initial_RP_decomp, creator_tag)
        self.patterns.append(null_pattern)

        # Create and store null RP structure...
        null_RP_structure = RP_structure(
            self.initial_RP_decomp,
            creator_tag,
            None,
            null_pattern,
        )
        null_pattern.RP_structures.append(null_RP_structure)
        self.RP_structures.append(null_RP_structure)

        # For each word in word list...
        for word in self.initial_RP_decomp.word_list:

            # Create and store a root.
            root_ptr = self.retrieve_or_create_root(word, creator_tag)

            # Associate it to null RP structure.
            if root_ptr not in null_RP_structure.roots:
                null_RP_structure.roots.append(root_ptr)
                root_ptr.RP_structures.append(null_RP_structure)

            # Add corresponding word analysis to data description...
            self.data_description.append(
                Word_RP_analysis(
                    self.initial_RP_decomp,
                    null_RP_structure,
                    root_ptr,
                    null_pattern,
                )
            )

            # Increment counts.
            self.increment(null_RP_structure, root_ptr, null_pattern)

            # Compute and store description length.
            self.null_model_length = self.description_length()

    def bootstrap_heuristics(
        self,
        min_root_length=1,
        min_vowels_in_pattern=1,
        min_num_roots=1,
        min_num_patterns=1,
    ):
        """Bootstrap RP analysis based on initial RP decomp"""
        creator_tag = 'BS'
        print "Bootstrapping RP analysis..."

        tmp_patterns_for_root = dict()
        tmp_roots_for_patterns = dict()

        # For each word
        for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
            word = self.initial_RP_decomp.word_list[word_idx]

            # Decompose it...
            word_root, word_pattern = \
                self.initial_RP_decomp.word_decomp[word_idx]

            # If candidate root and pattern are long enough...
            if (
                len(word_root) >= min_root_length and
                len(word_pattern) >= min_vowels_in_pattern + len(word_root)
            ):

                # Store tmp root and pattern, associate them...
                try:
                    tmp_patterns_for_root[word_root].append(word_pattern)
                except KeyError:
                    tmp_patterns_for_root[word_root] = [word_pattern]

        # For each tmp root...
        for root in tmp_patterns_for_root:

            # Get list of unique associated patterns...
            patterns = tuple(sorted(set(tmp_patterns_for_root[root])))

            # If pattern list size > 1...
            if len(patterns) >= min_num_patterns:

                # Associate root to pattern list...
                try:
                    tmp_roots_for_patterns[patterns].append(root)
                except KeyError:
                    tmp_roots_for_patterns[patterns] = [root]

        # For each pattern list
        for patterns in tmp_roots_for_patterns:

            # If number of associated roots > 1...
            if len(tmp_roots_for_patterns[patterns]) >= min_num_roots:

                # Create new RP structure, roots and patterns...
                new_RP_structure = RP_structure(
                    self.initial_RP_decomp,
                    creator_tag,
                )
                for root in tmp_roots_for_patterns[patterns]:
                    root_ptr = self.retrieve_or_create_root(root, creator_tag)
                    root_ptr.RP_structures.append(new_RP_structure)
                    new_RP_structure.roots.append(root_ptr)
                for pattern in patterns:
                    pattern_ptr = self.retrieve_or_create_pattern(
                        pattern,
                        creator_tag
                    )
                    pattern_ptr.RP_structures.append(new_RP_structure)
                    new_RP_structure.patterns.append(pattern_ptr)
                self.RP_structures.append(new_RP_structure)

        # For each word
        for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
            word = self.initial_RP_decomp.word_list[word_idx]

            # Decompose it, retrieve root and pattern...
            word_root, word_pattern = \
                self.initial_RP_decomp.word_decomp[word_idx]
            new_root = self.retrieve_root(word_root)
            new_pattern = self.retrieve_pattern(word_pattern)

            # Skip if no analysis or pattern is empty...
            if new_root is None or new_pattern is None:
                continue

            # Get and decrement old root, pattern and RP structure...
            self.decrement(
                self.data_description[word_idx].RP_structure,
                self.data_description[word_idx].root,
                self.data_description[word_idx].pattern,
            )

            # Update word analysis...
            self.data_description[word_idx].RP_structure   \
                = new_root.RP_structures[0]
            self.data_description[word_idx].root = new_root
            self.data_description[word_idx].pattern = new_pattern

            # Increment count of new root, pattern, RP structure...
            self.increment(
                new_root.RP_structures[0],
                new_root,
                new_pattern,
            )

    def extend_known_RP_structures(self, iter=1):
        """Find unanalyzed words that can be added to existing RP structures"""
        creator_tag = 'EKRPS' + str(iter)
        print "Extending known RP structures..."

        # For each RP structure (sorted by robustness)
        for RP_structure in sorted(
            self.RP_structures[1:], key=lambda r: r.robustness(), reverse=True
        ):

            # Cache current RP analysis.
            cached_RP_analysis = copy.deepcopy(self)

            # Initializations...
            tmp_patterns_for_root = dict()
            roots_to_be_added = list()

            # For each word
            for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
                word = self.initial_RP_decomp.word_list[word_idx]

                # Skip if analyzed.
                if self.data_description[word_idx].pattern.string is not None:
                    continue

                # Decompose word...
                word_root, word_pattern = \
                    self.initial_RP_decomp.word_decomp[word_idx]

                # Store tmp root and pattern, associate them...
                try:
                    tmp_patterns_for_root[word_root].append(word_pattern)
                except KeyError:
                    tmp_patterns_for_root[word_root] = [word_pattern]

            # For each tmp root...
            for root in tmp_patterns_for_root:

                # Skip if known.
                if self.retrieve_root(root) is not None:
                    continue

                # Get list of unique associated patterns...
                patterns = tuple(sorted(set(tmp_patterns_for_root[root])))

                # If it occurs with each pattern in this RP structure...
                root_can_be_added = True
                for pattern in RP_structure.patterns:
                    if pattern.string not in patterns:
                        root_can_be_added = False
                        break

                # Add it to list of roots to be added...
                if root_can_be_added == True:
                    roots_to_be_added.append(root)

            # Skip if no roots to be added...
            if len(roots_to_be_added) == 0:
                continue

            # Create roots and add them to RP structure...
            for root in roots_to_be_added:
                root_ptr = self.retrieve_or_create_root(root, creator_tag)
                root_ptr.RP_structures = [RP_structure]
                RP_structure.roots.append(root_ptr)

            # For each word
            for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
                word = self.initial_RP_decomp.word_list[word_idx]

                # Skip if analyzed.
                if self.data_description[word_idx].pattern.string is not None:
                    continue

                # Decompose it.
                word_root, word_pattern = \
                    self.initial_RP_decomp.word_decomp[word_idx]
                new_root = self.retrieve_root(word_root)
                new_pattern = self.retrieve_pattern(word_pattern)

                # If newly analyzed word...
                if (
                    word_root in roots_to_be_added and
                    new_pattern in RP_structure.patterns
                ):

                    # Decrement count of old root, pattern and RP structure...
                    self.decrement(
                        self.data_description[word_idx].RP_structure,
                        self.data_description[word_idx].root,
                        self.data_description[word_idx].pattern,
                    )

                    # Update word analysis...
                    self.data_description[word_idx].RP_structure   \
                        = new_root.RP_structures[0]
                    self.data_description[word_idx].root = new_root
                    self.data_description[word_idx].pattern = new_pattern

                    # Increment count of new root, pattern, RP structure...
                    self.increment(
                        new_root.RP_structures[0],
                        new_root,
                        new_pattern,
                    )

            # If description length is reduced, validate change...
            if cached_RP_analysis.description_length()  \
                > self.description_length():
                cached_RP_analysis = copy.deepcopy(self)
            else:
                self.revert(cached_RP_analysis)

    def extend_fragments_known_RP_structures(
        self,
        min_num_patterns=1,
        iter=1,
    ):
        """Find roots that are associated with a subset of the patterns of
        an existing RP structure
        """
        creator_tag = 'EFKRPS' + str(iter)
        print "Extending fragments of known RP structures..."

        # For each RP structure (sorted by robustness)
        for known_RP_structure in sorted(
            self.RP_structures[1:], key=lambda r: r.robustness(), reverse=True
        ):

            # Cache current RP analysis.
            cached_RP_analysis = copy.deepcopy(self)

            # Initializations...
            tmp_patterns_for_root = dict()
            roots_to_be_created = list()
            tmp_roots_for_patterns = dict()

            # For each word
            for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
                word = self.initial_RP_decomp.word_list[word_idx]

                # Skip if analyzed.
                if self.data_description[word_idx].pattern.string is not None:
                    continue

                # Decompose word...
                word_root, word_pattern = \
                    self.initial_RP_decomp.word_decomp[word_idx]

                # Store tmp root and pattern, associate them...
                try:
                    tmp_patterns_for_root[word_root].append(word_pattern)
                except KeyError:
                    tmp_patterns_for_root[word_root] = [word_pattern]

            # For each tmp root...
            for root in tmp_patterns_for_root:

                # Skip if known.
                if self.retrieve_root(root) is not None:
                    continue

                # Get list of unique associated patterns...
                patterns = tuple(sorted(set(tmp_patterns_for_root[root])))

                # Get list of patterns from this RP structure that occur with
                # this root...
                known_patterns = tuple(
                    sorted(
                        [
                            p1 for p1 in patterns
                            if p1 in [
                                p2.string for p2 in known_RP_structure.patterns
                            ]
                        ]
                    )
                )

                # If this list is large enough...
                if len(known_patterns) >= min_num_patterns:

                    # Add the root to list of roots to be created...
                    roots_to_be_created.append(root)

                    # Associate it to pattern list...
                    try:
                        tmp_roots_for_patterns[known_patterns].append(root)
                    except KeyError:
                        tmp_roots_for_patterns[known_patterns] = [root]

            # Skip if no roots to be created...
            if len(roots_to_be_created) == 0:
                continue

            # For each pattern list
            for patterns in tmp_roots_for_patterns:

                # Create new RP structure, roots and patterns...
                new_RP_structure = RP_structure(
                    self.initial_RP_decomp,
                    creator_tag,
                )
                for root in tmp_roots_for_patterns[patterns]:
                    root_ptr = self.retrieve_or_create_root(root, creator_tag)
                    root_ptr.RP_structures.append(new_RP_structure)
                    new_RP_structure.roots.append(root_ptr)
                for pattern in patterns:
                    pattern_ptr = self.retrieve_or_create_pattern(
                        pattern,
                        creator_tag
                    )
                    pattern_ptr.RP_structures.append(new_RP_structure)
                    new_RP_structure.patterns.append(pattern_ptr)
                self.RP_structures.append(new_RP_structure)

                # For each word
                for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
                    word = self.initial_RP_decomp.word_list[word_idx]

                    # Skip if analyzed.
                    if self.data_description[word_idx].pattern.string is not None:
                        continue

                    # Decompose it.
                    word_root, word_pattern = \
                        self.initial_RP_decomp.word_decomp[word_idx]
                    new_root = self.retrieve_root(word_root)
                    pattern = self.retrieve_pattern(word_pattern)

                    # If newly analyzed word...
                    if (
                        new_root in new_RP_structure.roots and
                        pattern in new_RP_structure.patterns
                    ):

                        # Decrement count of old root, pattern and RP structure...
                        self.decrement(
                            self.data_description[word_idx].RP_structure,
                            self.data_description[word_idx].root,
                            self.data_description[word_idx].pattern,
                        )

                        # Update word analysis...
                        self.data_description[word_idx].RP_structure   \
                            = new_RP_structure
                        self.data_description[word_idx].root = new_root
                        self.data_description[word_idx].pattern = pattern

                        # Increment count of new root, pattern, RP structure...
                        self.increment(
                            new_RP_structure,
                            new_root,
                            pattern,
                        )

            # If description length is reduced, validate change...
            if cached_RP_analysis.description_length()  \
                > self.description_length():
                cached_RP_analysis = copy.deepcopy(self)
            else:
                self.revert(cached_RP_analysis)

    def extend_known_roots_to_known_patterns(self, iter=1):
        """Find words that are made up of known roots and known patterns
        """
        creator_tag = 'EKRKP' + str(iter)
        print "Extending known roots to known patterns..."

        # Initializations...
        tmp_patterns_for_root = dict()
        tmp_roots_for_patterns = dict()

        # For each word
        for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
            word = self.initial_RP_decomp.word_list[word_idx]

            # Skip if analyzed.
            if self.data_description[word_idx].pattern.string is not None:
                continue

            # Decompose word...
            word_root, word_pattern = \
                self.initial_RP_decomp.word_decomp[word_idx]
            root_ptr = self.retrieve_root(word_root)
            pattern_ptr = self.retrieve_pattern(word_pattern)

            # If root and pattern are known...
            if root_ptr and pattern_ptr:

                # Store tmp root and pattern, associate them...
                try:
                    tmp_patterns_for_root[root_ptr].append(pattern_ptr)
                except KeyError:
                    tmp_patterns_for_root[root_ptr] = [pattern_ptr]

        # For each tmp root...
        for root_ptr in tmp_patterns_for_root:

            # Get list of unique associated patterns (including those from
            # the original RP structure)...
            patterns = tuple(
                sorted(
                    set(tmp_patterns_for_root[root_ptr]).union(
                        set(root_ptr.RP_structures[0].patterns)
                    )
                )
            )

            # Associate root to pattern list...
            try:
                tmp_roots_for_patterns[patterns].append(root_ptr)
            except KeyError:
                tmp_roots_for_patterns[patterns] = [root_ptr]

        # For each pattern list
        for patterns in tmp_roots_for_patterns:

            # Cache current RP analysis.
            cached_RP_analysis = copy.deepcopy(self)

            # Create new RP structure, associate roots and patterns...
            new_RP_structure = RP_structure(
                self.initial_RP_decomp,
                creator_tag,
            )
            for root_ptr in tmp_roots_for_patterns[patterns]:
                root_ptr.RP_structures.append(new_RP_structure)
                new_RP_structure.roots.append(root_ptr)
            for pattern_ptr in patterns:
                pattern_ptr.RP_structures.append(new_RP_structure)
                new_RP_structure.patterns.append(pattern_ptr)
            self.RP_structures.append(new_RP_structure)

            # For each word
            for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
                word = self.initial_RP_decomp.word_list[word_idx]

                # Decompose it.
                word_root, word_pattern = \
                    self.initial_RP_decomp.word_decomp[word_idx]
                root_ptr = self.retrieve_root(word_root)
                pattern_ptr = self.retrieve_pattern(word_pattern)

                # If newly analyzed word...
                if (
                    root_ptr in new_RP_structure.roots and
                    pattern_ptr in new_RP_structure.patterns
                ):

                    # Decrement count of old root, pattern and RP structure...
                    self.decrement(
                        self.data_description[word_idx].RP_structure,
                        self.data_description[word_idx].root,
                        self.data_description[word_idx].pattern,
                    )

                    # Update word analysis...
                    self.data_description[word_idx].RP_structure   \
                        = new_RP_structure
                    self.data_description[word_idx].root = root_ptr
                    self.data_description[word_idx].pattern = pattern_ptr

                    # Increment count of new root, pattern, RP structure...
                    self.increment(
                        new_RP_structure,
                        root_ptr,
                        pattern_ptr,
                    )

            # If description length is reduced, validate change...
            if cached_RP_analysis.description_length()  \
                > self.description_length():
                cached_RP_analysis = copy.deepcopy(self)
            else:
                self.revert(cached_RP_analysis)

    def extend_known_roots(self, min_robustness, min_pattern_count=1, iter=1):
        """Find words that are made up of robust known roots"""
        creator_tag = 'EKR' + str(iter)
        print "Extending known roots..."

        # Initialization...
        tmp_patterns_for_root = dict()
        tmp_roots_for_patterns = dict()
        candidate_patterns = list()
        candidate_roots = list()

        # For each word
        for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
            word = self.initial_RP_decomp.word_list[word_idx]

            # Skip if analyzed.
            if self.data_description[word_idx].pattern.string is not None:
                continue

            # Decompose word...
            word_root, word_pattern = \
                self.initial_RP_decomp.word_decomp[word_idx]

            # Store tmp root and pattern, associate them...
            try:
                tmp_patterns_for_root[word_root].append(word_pattern)
            except KeyError:
                tmp_patterns_for_root[word_root] = [word_pattern]

        # For each root...
        for known_root in self.roots:

            # Skip if never associated with unanalyzed words...
            if known_root.string not in tmp_patterns_for_root:
                continue

            # Skip if always associated with unanalyzed words...
            if known_root.creator == 'INIT':
                continue

            # If the robustness of this root's RP structure is high enough...
            if known_root.RP_structures[0].robustness() >= min_robustness:

                # Add the root and patterns found with this root in unanalyzed
                # words to lists of candidates..
                candidate_patterns.extend(
                    tmp_patterns_for_root[known_root.string]
                )
                candidate_roots.append(known_root.string)

        # Keep only candidate patterns that are frequent enough...
        candidate_patterns = [
            p for p in candidate_patterns
            if candidate_patterns.count(p) >= min_pattern_count
        ]

        # From patterns associated with candidate roots, keep only candidates,
        candidate_roots_to_remove = list()
        for root in candidate_roots:
            tmp_patterns_for_root[root] = [
                p for p in tmp_patterns_for_root[root]
                if p in candidate_patterns
            ]
            if tmp_patterns_for_root[root] == []:
                candidate_roots_to_remove.append(root)

        # Remove candidate roots that have no candidate patterns associated...
        candidate_roots = [
            r for r in candidate_roots
            if r not in candidate_roots_to_remove
        ]

        # For each candidate root...
        for root in candidate_roots:

            # Get list of unique associated patterns (including those from
            # the original RP structure)...
            root_ptr = self.retrieve_root(root)
            patterns = tuple(
                sorted(
                    set(tmp_patterns_for_root[root]).union(
                        set(
                            p.string
                            for p in root_ptr.RP_structures[0].patterns
                        )
                    )
                )
            )

            # Associate root to pattern list...
            try:
                tmp_roots_for_patterns[patterns].append(root)
            except KeyError:
                tmp_roots_for_patterns[patterns] = [root]

        # For each pattern list
        for patterns in tmp_roots_for_patterns:

            # Cache current RP analysis.
            cached_RP_analysis = copy.deepcopy(self)

            # Create new RP structure, associate roots and patterns...
            new_RP_structure = RP_structure(
                self.initial_RP_decomp,
                creator_tag,
            )
            for root in tmp_roots_for_patterns[patterns]:
                root_ptr = self.retrieve_root(root)
                root_ptr.RP_structures.append(new_RP_structure)
                new_RP_structure.roots.append(root_ptr)
            for pattern in patterns:
                pattern_ptr = self.retrieve_or_create_pattern(
                    pattern,
                    creator_tag
                )
                pattern_ptr.RP_structures.append(new_RP_structure)
                new_RP_structure.patterns.append(pattern_ptr)
            self.RP_structures.append(new_RP_structure)

            # For each word
            for word_idx in xrange(len(self.initial_RP_decomp.word_list)):
                word = self.initial_RP_decomp.word_list[word_idx]

                # Decompose it.
                word_root, word_pattern = \
                    self.initial_RP_decomp.word_decomp[word_idx]
                root_ptr = self.retrieve_root(word_root)
                pattern_ptr = self.retrieve_pattern(word_pattern)

                # If newly analyzed word...
                if (
                    root_ptr in new_RP_structure.roots and
                    pattern_ptr in new_RP_structure.patterns
                ):

                    # Decrement count of old root, pattern and RP structure...
                    self.decrement(
                        self.data_description[word_idx].RP_structure,
                        self.data_description[word_idx].root,
                        self.data_description[word_idx].pattern,
                    )

                    # Update word analysis...
                    self.data_description[word_idx].RP_structure   \
                        = new_RP_structure
                    self.data_description[word_idx].root = root_ptr
                    self.data_description[word_idx].pattern = pattern_ptr

                    # Increment count of new root, pattern, RP structure...
                    self.increment(
                        new_RP_structure,
                        root_ptr,
                        pattern_ptr,
                    )

            # If description length is reduced, validate change...
            if cached_RP_analysis.description_length()  \
                > self.description_length():
                cached_RP_analysis = copy.deepcopy(self)
            else:
                self.revert(cached_RP_analysis)

    def revert(self, RP_analysis):
        self.initial_RP_decomp = RP_analysis.initial_RP_decomp
        self.roots = RP_analysis.roots
        self.patterns = RP_analysis.patterns
        self.RP_structures = RP_analysis.RP_structures
        self.data_description = RP_analysis.data_description
        self.null_model_length = RP_analysis.null_model_length

    def increment(self, RP_structure, root, pattern):
        """Increment count of RP structure, root and pattern"""
        root.count += 1
        pattern.count += 1
        RP_structure.count += 1
        RP_structure.root_count[root]   \
            = RP_structure.root_count.get(root, 0) + 1
        RP_structure.pattern_count[pattern]     \
            = RP_structure.pattern_count.get(pattern, 0) + 1

    def decrement(self, RP_structure, root, pattern):
        """Decrement count of RP structure, root and pattern"""
        RP_structure.root_count[root] -= 1
        if RP_structure.root_count[root] == 0:
            del RP_structure.root_count[root]
            RP_structure.roots.remove(root)
            root.RP_structures.remove(RP_structure)
        root.count -= 1
        if root.count == 0:
            self.roots.remove(root)
        RP_structure.pattern_count[pattern] -= 1
        if RP_structure.pattern_count[pattern] == 0:
            del RP_structure.pattern_count[pattern]
            RP_structure.patterns.remove(pattern)
            pattern.RP_structures.remove(RP_structure)
        pattern.count -= 1
        if pattern.count == 0 and pattern.creator != 'INIT':
            self.patterns.remove(pattern)
        RP_structure.count -= 1
        if RP_structure.count == 0 and RP_structure.creator != 'INIT':
            self.RP_structures.remove(RP_structure)

    def description_length(self):
        """Compute overall description length of RP analysis"""
        return (
              self.description_length_model()
            + self.description_length_data()
        )

    def description_length_model(self):
        """Compute description length of model"""
        return (
              self.description_length_roots()
            + self.description_length_patterns()
            + self.description_length_RP_structures()
        )

    def description_length_roots(self):
        """Compute description length of list of roots"""
        total_description_length = 0
        try:
            total_description_length += math.log(len(self.roots), 2)
        except ValueError:
            pass
        for root in self.roots:
            total_description_length += root.description_length()
        return total_description_length

    def description_length_patterns(self):
        """Compute description length of list of patterns"""
        total_description_length = 0
        try:
            total_description_length += math.log(len(self.patterns), 2)
        except ValueError:
            pass
        for pattern in self.patterns:
            total_description_length += pattern.description_length()
        return total_description_length

    def description_length_RP_structures(self):
        """Compute description length of list of RP_structures"""
        total_description_length = 0
        try:
            total_description_length += math.log(len(self.RP_structures), 2)
        except ValueError:
            pass
        for RP_structure in self.RP_structures:
            total_description_length += RP_structure.description_length()
        return total_description_length

    def description_length_data(self):
        """Compute description length of data description"""
        total_description_length = 0
        try:
            total_description_length += list_description_length(
                len(self.initial_RP_decomp.word_list)
            )
        except ValueError:
            pass
        for word_analysis in self.data_description:
            total_description_length += word_analysis.description_length()
        return total_description_length

    def display(self):
        """Display the details of a RP analysis"""
        DL_roots            = self.description_length_roots()
        DL_patterns         = self.description_length_patterns()
        DL_RP_structures    = self.description_length_RP_structures()
        DL_data_description = self.description_length_data()
        total_DL = DL_roots+DL_patterns+DL_RP_structures+DL_data_description
        delta = self.null_model_length - total_DL
        relative_delta = delta / self.null_model_length
        output  = ''
        output += '-' * 78 + '\n* DESCRIPTION LENGTH ANALYSIS\n' + '-'*78 + '\n'
        output += 'List of roots:\t%.2f bits\n' % DL_roots
        output += 'List of patterns:\t%.2f bits\n' % DL_patterns
        output += 'List of RP structures:\t%.2f bits\n' % DL_RP_structures
        output += 'Compressed corpus:\t%.2f bits\n' % DL_data_description
        output += 'Total:\t%.2f bits\n' % total_DL
        output += 'Delta:\t%.2f bits\n' % delta
        output += 'Relative delta:\t%.2f%%\n' % (relative_delta * 100)
        output += '\n'
        output += '-' * 78 + '\n* COMPRESSED CORPUS\n' + '-' * 78 + '\n'
        output += 'WORD\tRP STRUCTURE\tDL\n'
        output += '\n'.join(
            a.display() for a in sorted(
                self.data_description, key=lambda a: a.compose()
            )
        ) + '\n'
        output += '\n'
        output += '-' * 78 + '\n* PHONOLOGICAL CATEGORIES\n' + '-' * 78 + '\n'
        output += 'C = {%s}\n' % ','.join(
            sorted(self.initial_RP_decomp.consonants)
        )
        output += 'V = {%s}\n' % ','.join(
            sorted(self.initial_RP_decomp.vowels)
        )
        output += '\n'
        output += '-' * 78 + '\n* ROOTS\n' + '-' * 78 + '\n'
        output += 'ROOT\tCOUNT\tRP STRUCTURES\tDL\tCREATOR\n'
        output += '\n'.join(
            r.display() for r in sorted(self.roots, key=lambda r: r.string)
        ) + '\n'
        output += '\n'
        output += '-' * 78 + '\n* PATTERNS\n' + '-' * 78 + '\n'
        output += 'PATTERN\tCOUNT\tRP STRUCTURES\tDL\tCREATOR\n'
        output += '\n'.join(
            p.display() for p in sorted(self.patterns, key=lambda p: p.string)
        ) + '\n'
        output += '\n'
        output += '-' * 78 + '\n* RP STRUCTURES\n' + '-' * 78 + '\n'
        output += ('-' * 78 + '\n').join(
            r.display() for r in sorted(
                self.RP_structures, key=lambda r: r.robustness(), reverse=True
            )
        ) + '\n'
        return output

    def retrieve_root(self, root):
        """Look for root and return pointer to it or None"""
        return next(
            (r for r in self.roots if r.string == root),
            None
        )

    def retrieve_or_create_root(self, root, creator_tag, init=0):
        """Look for root and return pointer to it or to a new one"""
        root_ptr = self.retrieve_root(root)
        if root_ptr is None:
            root_ptr = Root(self.initial_RP_decomp, creator_tag, root)
            self.roots.append(root_ptr)
            root_ptr.count = init
        return root_ptr

    def retrieve_pattern(self, pattern):
        """Look for pattern and return pointer to it or None"""
        return next(
            (p for p in self.patterns if p.string == pattern),
            None
        )

    def retrieve_or_create_pattern(self, pattern, creator_tag, init=0):
        """Look for pattern and return pointer to it or to a new one"""
        pattern_ptr = self.retrieve_pattern(pattern)
        if pattern_ptr is None:
            pattern_ptr = Pattern(self.initial_RP_decomp, creator_tag, pattern)
            self.patterns.append(pattern_ptr)
            pattern_ptr.count = init
        return pattern_ptr

    def retrieve_RP_structure(self, root_ptr, pattern_ptr):
        """Look for RP structure and return pointer to it or None"""
        return next(
            (
                s for s in self.RP_structures if (
                        root_ptr    in s.roots
                    and pattern_ptr in s.patterns
                )
            ),
            None
        )

    def retrieve_or_create_RP_structure(
        self,
        root_ptr,
        pattern_ptr,
        creator_tag,
        init=0
    ):
        """Look for RP structure and return pointer to it or to a new one"""
        RP_structure_ptr = self.retrieve_RP_structure(
            root_ptr,
            pattern_ptr,
        )
        if RP_structure_ptr is None:
            RP_structure_ptr = RP_structure(
                self.initial_RP_decomp,
                creator_tag,
                root_ptr,
                pattern_ptr,
            )
            self.RP_structures.append(RP_structure_ptr)
        RP_structure_ptr.count = init
        return RP_structure_ptr




class Morpheme(object):
    """Abstraction for root or pattern"""
    def __init__(self, string=None):
        self.string = string
        self.count  = 0

    def __str__(self):
        return self.string or 'NULL'

    def display(self):
        output  = '%s\t%i\t{' % (self, self.count)
        output += ','.join(str(r.id) for r in self.RP_structures)
        output += '}\t%.2f\t' % self.description_length()
        output += self.creator
        return output


class Root(Morpheme):
    """Class for representing a root"""
    class_counter = 0
    def __init__(
        self,
        initial_RP_decomp,
        creator,
        string=None,
        RP_structures=None
    ):
        super(Root, self).__init__(string)
        self.initial_RP_decomp  = initial_RP_decomp
        self.creator            = creator
        self.RP_structures      = RP_structures or list()
        self.id                 = Root.class_counter
        Root.class_counter += 1

    def description_length(self):
        alphabet_size = len(self.initial_RP_decomp.consonants) \
                      + len(self.initial_RP_decomp.vowels)
        return (
              len(self.string)
            * math.log(alphabet_size, 2)
            + list_description_length(len(self.string))
        )


class Pattern(Morpheme):
    """Class for representing a pattern"""
    class_counter = 0
    def __init__(
        self,
        initial_RP_decomp,
        creator,
        string=None,
        RP_structures=None
    ):
        super(Pattern, self).__init__(string)
        self.initial_RP_decomp  = initial_RP_decomp
        self.creator            = creator
        self.RP_structures      = RP_structures or list()
        self.id                 = Pattern.class_counter
        Pattern.class_counter += 1
    def description_length(self):
        alphabet_size = len(self.initial_RP_decomp.vowels) \
                      + 1
        try:
            return (
                  len(self.string)
                * math.log(alphabet_size, 2)
                + list_description_length(len(self.string))
            )
        except TypeError:
            return 0


class RP_structure(object):
    """Class for representing an RP structure"""
    class_counter = 0

    def __init__(
        self,
        initial_RP_decomp,
        creator,
        root_ptr=None,
        pattern_ptr=None
    ):
        self.initial_RP_decomp  = initial_RP_decomp
        self.creator            = creator
        self.roots              = list()
        self.root_count         = dict()
        self.patterns           = list()
        self.pattern_count      = dict()
        if root_ptr:
            self.roots.append(root_ptr)
        if pattern_ptr:
            self.patterns.append(pattern_ptr)
        self.count  = 0
        self.id = RP_structure.class_counter
        RP_structure.class_counter += 1

    def robustness(self):
        sum_len_roots    = sum(len(r.string) for r in self.roots)
        sum_len_patterns = sum(len(p.string) for p in self.patterns if p.string)
        if sum_len_patterns > 0:
            return (len(self.roots) - 1) * sum_len_patterns - sum_len_roots
        else:
            return None    # Default for NULL RP structure (NA)

    def description_length(self):
        total = 0
        try:
            total += list_description_length(len(self.roots))
        except ValueError:
            pass
        try:
            total += list_description_length(len(self.patterns))
        except ValueError:
            pass
        num_words = len(self.initial_RP_decomp.word_list)
        for root in self.roots:
            prob = root.count / num_words
            try:
                total -= math.log(prob, 2)
            except ValueError:
                pass
        for pattern in self.patterns:
            prob = pattern.count / num_words
            try:
                total -= math.log(prob, 2)
            except ValueError:
                pass
        return total

    def display(self):
        output  = 'ID\tROBUSTNESS\tCOUNT\tDL\tCREATOR\n'
        output += '%i\t%s\t%i\t%.2f\t%s\n' % (
            self.id,
            self.robustness() or "NA",
            self.count,
            self.description_length(),
            self.creator,
        )
        output += '\t' + '-' * 20 + '\n\tROOT\tCOUNT\tDL\n'
        num_words = len(self.initial_RP_decomp.word_list)
        for root in sorted(
            self.roots,
            key=lambda r: self.root_count[r],
            reverse=True,
        ):
            prob = root.count / num_words
            try:
                description_length = -math.log(prob, 2)
            except ValueError:
                description_length = 0
            output += '\t%s\t%i\t%.2f\n' % (
                root,
                self.root_count[root],
                description_length,
            )
        output += '\t' + '-' * 20 + '\n\tPATTERN\tCOUNT\tDL\n'
        for pattern in sorted(
            self.patterns,
            key=lambda p: self.pattern_count[p],
            reverse=True,
        ):
            prob = pattern.count / num_words
            try:
                description_length = -math.log(prob, 2)
            except ValueError:
                description_length = 0
            if description_length == -0:
                description_length = 0
            output += '\t%s\t%i\t%.2f\n' % (
                pattern,
                self.pattern_count[pattern],
                description_length,
            )
        return output


class Word_RP_analysis(object):
    """Class for representing the RP analysis of a word"""
    def __init__(
        self,
        initial_RP_decomp,
        RP_structure_ptr=None,
        root_ptr=None,
        pattern_ptr=None
    ):
        self.initial_RP_decomp  = initial_RP_decomp
        self.RP_structure       = RP_structure_ptr
        self.root               = root_ptr
        self.pattern            = pattern_ptr

    def compose(self):
        try:
            return self.pattern.string.replace(ins_symbol, '%c') % tuple(
                self.root.string
            )
        except AttributeError:
            return self.root.string

    def description_length(self):
        total = 0
        num_words = len(self.initial_RP_decomp.word_list)
        total -= math.log(self.RP_structure.count / num_words, 2)
        total -= math.log(
            (
                self.RP_structure.root_count[self.root]
              / self.RP_structure.count
            ),
            2
        )
        total -= math.log(
            (
                self.RP_structure.pattern_count[self.pattern]
              / self.RP_structure.count
            ),
            2
        )
        return total

    def display(self):
        if self.pattern.string is not None:
            word = self.compose()
        else:
            word = self.root
        return '%s\t%i\t%.2f' % (
            word,
            self.RP_structure.id,
            self.description_length(),
        )

def list_description_length(num_items):
    return math.log(num_items, 2)


if __name__ == "__main__":
    main()

