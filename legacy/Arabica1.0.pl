#!/usr/bin/perl -w

use strict;

################################################################################
#
# FILE: Arabica.pl - Copyright 2007 Aris Xanthos.
#
# This file is part of Arabica 1.0
#
# Arabica 1.0 is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------
#
# SYNOPSIS: This is a program for the unsupervised analysis of root-and-pattern
# morphology. It analyzes a set of words in terms of 'RP structures', using
# Sukhotin's algorithm for vowel identification. It is inspired by the MDL
# approach to morphological learning developed by John Goldsmith (2001, 2006)
# and implemented in the software Linguistica (http://linguistica.uchicago.edu).
#
#-------------------------------------------------------------------------------
#
# INPUT: - A text file with 1 word per line (words can be repeated or not)
#        - The file Config_Arabica.txt where parameters are set
#
#-------------------------------------------------------------------------------
#
# OUTPUT: A text file summarizing the morphology that was learned
#
################################################################################





################################################################################
# GLOBAL VARIABLES DEFINITIONS

my( %Parameters, %Corpus, @RPStructureRobustness, @RPStructureFrequency );
my( %Consonants, %Vowels, %ListOfRoots, %ListOfPatterns, @ListOfRPStructures );
my( %WordsToRPStructures, %RootsToRPStructures, %PatternsToRPStructures );
my( $NumWords, %WordsToPatterns, $UnanalyzedWordTokens, $UnanalyzedWordTypes );
my( %WordsToRoots, %RootCreator, %RPDescriptionLength );

################################################################################





################################################################################
# MAIN PROGRAM:



#-------------------------------------------------------------------------------
# Read input files:

ReadParameters() || die "Parameters file error";

$NumWords = ReadCorpus() || die "Input file error";

#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# RP analysis:

# Infer consonants and vowels
Sukhotin();

# Check for notation conflicts
CheckNotation() || die "Notation conflict: please select another symbol for insertion slots.\n";

# Build initial RP structuress
RPBootstrapping();

# Extend known RP structures
if( $Parameters{ExtendKnownRPStructures} eq 'true' ) {
  ExtendKnownRPStructures();
}

# Compute description length of RP analysis
RPDescriptionLength();

#-------------------------------------------------------------------------------



OutputResults() || die "Output file error";



# END OF MAIN PROGRAM.
################################################################################





#-------------------------------------------------------------------------------
# FUNCTION: ReadParameters
#-------------------------------------------------------------------------------
# SYNOPSIS: Reads parameters from file Config_Arabica.txt and stores them in
# the %Parameters hash.
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: 1 if the parameters are read correctly, 0 otherwise
#-------------------------------------------------------------------------------

sub ReadParameters {

  my( $Line );

  open( INPUT, "Config_Arabica.txt" ) || return 0;

  # Set default parameters
  $Parameters{SamplingRate} = 1;
  $Parameters{MinimumNumberOfRoots} = 2;
  $Parameters{MinimumNumberOfPatterns} = 2;
  $Parameters{MinimumNumberOfPhonemesPerRoot} = 2;
  $Parameters{MinimumNumberOfPhonemesPerPattern} = 1;
  $Parameters{PhonemeCategoryInPatterns} = 'vowels';
  $Parameters{ExtendKnownRPStructures} = 'true';
  $Parameters{InsertionSlotSymbol} = '-';

  while( $Line = <INPUT> ) {
     chomp $Line;

    # Leave aside comments (lines that begin with '#') and empty lines
    unless( $Line =~ /^#/ || $Line eq "" ) {

      # If line has correct syntax (see Config_Arabica.txt)
      if( $Line =~ /^([^= ]+) *= *([^= ]+)$/ ) {
        $Parameters{$1} = $2; # Set value in %Parameters hash
      }
      else { return 0; }

    }
  }
  close INPUT;

  # Check that sampling rate is > 0 and <= 1
  if( $Parameters{SamplingRate} <= 0 || $Parameters{SamplingRate} > 1 ) {
    print "Parameter SamplingRate is out of range ]0,1] and will be set to 1\n";
    $Parameters{SamplingRate} = 1;
  }

  # Check that the symbol for insertion slots is not empty
  if( $Parameters{InsertionSlotSymbol} eq "" ) {
    print "Parameter InsertionSlotSymbol was not specified and will be set to '-'\n";
    $Parameters{InsertionSlotSymbol} = '-';
  }

  # Check that input file name is specified before exit
  if( defined $Parameters{InputFileName} ) { return 1; }
  else { return 0; }
}



#-------------------------------------------------------------------------------
# FUNCTION: ReadCorpus
#-------------------------------------------------------------------------------
# SYNOPSIS: Reads the specified corpus according to parameters and stores
# words in the %Corpus hash with their frequency.
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: The number of words that were read.
#-------------------------------------------------------------------------------

sub ReadCorpus {

  my( $Line, $CorpusRead );

  $CorpusRead = 0;

  open( INPUT, $Parameters{InputFileName} ) || return 0;

  while( $Line = <INPUT> ) {
    chomp $Line;

    # If the line is a string of characters without separators
    if( $Line =~ /^\S+$/ ){

      # If it is sampled
      if( rand() <= $Parameters{SamplingRate} ) {

        $Corpus{$Line}++;

        $CorpusRead++;

      }

    }

  }

  close INPUT;

  return $CorpusRead;
}



#-------------------------------------------------------------------------------
# FUNCTION: CheckNotation
#-------------------------------------------------------------------------------
# SYNOPSIS: Ensures that the symbol selected for insertion slots does not occur
# in the data.
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: 1 if there is no conflict, 0 otherwise.
#-------------------------------------------------------------------------------

sub CheckNotation {

  if( exists $Consonants{$Parameters{InsertionSlotSymbol}} ||
      exists $Vowels{$Parameters{InsertionSlotSymbol}} ) { return 0; }

  else{ return 1; }
}



#-------------------------------------------------------------------------------
# FUNCTION: RPBootstrapping
#-------------------------------------------------------------------------------
# SYNOPSIS: Analyzes each word into the conjunction of a consonantal root and
# a pattern (insertion slots and vowels). Stores each such component into the
# corresponding list, and creates an RP structure, i.e. a pair of lists of
# pointers to roots and patterns.
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: none
#-------------------------------------------------------------------------------

sub RPBootstrapping {

  my( $Word, %Analysis, $RPStructureRef, $NumberOfRPStructures, $RPStructureIndex );
  my( $RPStructureIndex2, @Array1, @Array2, $Root, $Pattern );

  $UnanalyzedWordTokens = 0;
  $UnanalyzedWordTypes = 0;

  $NumberOfRPStructures = 0;

  # For each word
  foreach $Word ( keys %Corpus ) {

    # Analysis into root-pattern components (the result is stored into
    # the %Analysis hash).
    Decompose( $Word, \%Analysis );

    # If the resulting root and pattern have the required number of phonemes

    if( length( $Analysis{root} ) >= $Parameters{MinimumNumberOfPhonemesPerRoot} &&
        NumberOfPhonemesInPattern( $Analysis{pattern} ) >= $Parameters{MinimumNumberOfPhonemesPerPattern} )  {

      # Store the components in the corresponding lists and increment their counts.
      $ListOfRoots{$Analysis{root}} += $Corpus{$Word};
      $ListOfPatterns{$Analysis{pattern}} += $Corpus{$Word};

      # Indicate at which stage this root was created (BS = bootstrapping)
      $RootCreator{$Analysis{root}} = 'BS';

      # Create a new RP structure with these components, and store a
      # reference to it in the list of RP structures.

      $RPStructureRef = NewRPStructure( [$Analysis{root}], [$Analysis{pattern}] );
      push( @ListOfRPStructures, $RPStructureRef );

      $RPStructureIndex = $NumberOfRPStructures++;

      # Assign to the word, root, and pattern a reference to the RP structure.
      $WordsToRPStructures{$Word} = $RPStructureIndex;
      AddValueToArray( \@{ $RootsToRPStructures{$Analysis{root}} }, $RPStructureIndex );
      AddValueToArray( \@{ $PatternsToRPStructures{$Analysis{pattern}} }, $RPStructureIndex );

    }

    # If the resulting root and pattern do NOT have the required number of phonemes
    else {

      # Assign to the word a reference to the null RP structure
      $WordsToRPStructures{$Word} = -1;

      # Increment the type and token counts of words in null RP structure.
      $UnanalyzedWordTokens += $Corpus{$Word};
      $UnanalyzedWordTypes++;

    }

  }

  # Collapse RP structures with identical roots

  # For each RP structure
  for( $RPStructureIndex = 0;
       $RPStructureIndex < $NumberOfRPStructures-1; $RPStructureIndex++ ) {

    @Array1 = keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} };

    # Loop over each other RP structure
    for( $RPStructureIndex2 = $RPStructureIndex+1;
         $RPStructureIndex2 < $NumberOfRPStructures; $RPStructureIndex2++ ) {

      @Array2 = keys %{ ${ $ListOfRPStructures[$RPStructureIndex2] }{roots} };

      # If they have identical roots (TO DO: relax identity requirement)
      if( CompareArrays( \@Array1, \@Array2 ) ) {

        # Collapse them
        CollapseRPStructures( $RPStructureIndex, $RPStructureIndex2 );

        # Delete the 2nd RP structure.
        DeleteRPStructure( $RPStructureIndex2 );

        $RPStructureIndex2--;
        $NumberOfRPStructures--;
      }

    }

  }

  # Collapse RP structures with identical patterns

  # For each RP structure
  for( $RPStructureIndex = 0;
       $RPStructureIndex < $NumberOfRPStructures-1; $RPStructureIndex++ ) {

    @Array1 = keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} };

    # Loop over each other RPStructure
    for( $RPStructureIndex2 = $RPStructureIndex+1;
         $RPStructureIndex2 < $NumberOfRPStructures; $RPStructureIndex2++ ) {

      @Array2 = keys %{ ${ $ListOfRPStructures[$RPStructureIndex2] }{patterns} };

      # If they have identical patterns (TO DO: relax identity requirement)
      if( CompareArrays( \@Array1, \@Array2 ) ) {

        # Collapse them
        CollapseRPStructures( $RPStructureIndex, $RPStructureIndex2 );

        # Delete the 2nd RP structure.
        DeleteRPStructure( $RPStructureIndex2 );

        $RPStructureIndex2--;
        $NumberOfRPStructures--;
      }

    }

  }

  # Keep only those RP structures that have at least the number of roots and
  # patterns specified in the parameters.

  # For each RP structure
  for( $RPStructureIndex = 0;
       $RPStructureIndex < $NumberOfRPStructures; $RPStructureIndex++ ) {

    # Delete it if it does not have the required number of roots and patterns.
    if( keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} } <
          $Parameters{MinimumNumberOfRoots} ||
        keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} } <
          $Parameters{MinimumNumberOfPatterns} ) {

      # Increment the type and token counts of words in null RP structure.
      $UnanalyzedWordTokens += RPStructureFrequency( $RPStructureIndex );
      $UnanalyzedWordTypes += RPStructureFrequency( $RPStructureIndex, 'types' );

      DeleteRPStructure( $RPStructureIndex );
      $RPStructureIndex--;
      $NumberOfRPStructures--;

    }
    # Else calculate its robustness and frequency and link words to roots and patterns.
    else {
      $RPStructureRobustness[$RPStructureIndex] = RPStructureRobustness( $RPStructureIndex );
      $RPStructureFrequency[$RPStructureIndex] = RPStructureFrequency( $RPStructureIndex );
      LinkWordsToRootsAndPatterns( $RPStructureIndex );
    }

  }

}





#-------------------------------------------------------------------------------
# FUNCTION: ExtendKnownRPStructures
#-------------------------------------------------------------------------------
# SYNOPSIS: This heuristic tries to find sets of unanalyzed words that fit into
# one of the RP structures that were previously discovered. Such sets are
# integrated in the corresponding RP structure if they result in a decrease
# in description length.
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub ExtendKnownRPStructures {

  my( $Word, $Root, $Pattern, @RobustnessOrder, $RobustnessIndex, %Analysis );
  my( $RPStructureIndex, $NumberOfRPStructures, $MatchesRPStructure );
  my( $MatchesPattern, $ComposedWord, %TempHashOfRoots, $Count );

  $NumberOfRPStructures = int (@ListOfRPStructures);

  # Construct an array with the order of robustness of RP structures.
  RPStructureRobustnessOrder( \@RobustnessOrder );

  # For each RP structure (sorted by robustness)
  for( $RPStructureIndex = 0; $RPStructureIndex < $NumberOfRPStructures; $RPStructureIndex++ ) {

    $RobustnessIndex = $RobustnessOrder[$RPStructureIndex];

    undef %TempHashOfRoots;

    # Scan each word in the corpus
    foreach $Word ( keys %Corpus ) {

      # If the word is unanalyzed, determine if it contains one of the
      # patterns of the RP structure AND an unknown root.
      if( $WordsToRPStructures{$Word} == -1) {

        Decompose( $Word, \%Analysis );

        if( defined $RootsToRPStructures{$Analysis{root}} == 1 ||
            defined $TempHashOfRoots{$Analysis{root}} == 1 ) { next; }

        $MatchesPattern = 0;

        foreach $Pattern
          (keys %{ ${ $ListOfRPStructures[$RobustnessIndex] }{patterns} }) {

          if( $Analysis{pattern} eq $Pattern ) {
            $MatchesPattern = 1;
            last;
          }
        }

        # If it does, check whether the corresponding root occurs with all the
        # patterns in this RP structure (in unanalyzed words).
        if( $MatchesPattern == 1 ) {

          $MatchesRPStructure = 1;

          $Count = 0;

          foreach $Pattern
            (keys %{ ${ $ListOfRPStructures[$RobustnessIndex] }{patterns} }) {

            $ComposedWord = Compose( $Analysis{root}, $Pattern );

            if( exists $Corpus{ $ComposedWord } == 0 ||
                $WordsToRPStructures{ $ComposedWord } != -1 )
            {
              $MatchesRPStructure = 0;
              last;
            }
            else { $Count += $Corpus{ $ComposedWord }; }
          }

          # If it does, store this root in a temporary hash
          if( $MatchesRPStructure == 1 ) { $TempHashOfRoots{$Analysis{root}} = $Count;}

        }
      }
    }

    # Evaluate the difference in description length resulting from the addition
    # of the collected roots to the RP structure we are examining. If it is
    # negative (i.e. if it saves bits), accept the modification.

    if( int keys %TempHashOfRoots > 0 &&
        DeltaDescriptionLength( $RobustnessIndex, \%TempHashOfRoots ) < 0 ) {

      foreach $Root (keys %TempHashOfRoots) {

        foreach $Pattern
          ( keys %{ ${ $ListOfRPStructures[$RobustnessIndex] }{patterns} } )
        {

          $ComposedWord = Compose( $Root, $Pattern );
          $Count = $Corpus{ $ComposedWord };

          # Increment the count of roots and patterns in RP structure
          ${ $ListOfRPStructures[$RobustnessIndex] }{roots}{$Root} += $Count;
          ${ $ListOfRPStructures[$RobustnessIndex] }{patterns}{$Pattern} += $Count;

          # Store the components in the corresponding lists and increment their counts.
          $ListOfRoots{$Root} += $Count;
          $ListOfPatterns{$Pattern} += $Count;

          # Indicate at which stage this root was created (EXT = extend known
          # RP structures)
          $RootCreator{$Root} = 'EXT';

          # Decrement the count of unanalyzed word types and tokens
          $UnanalyzedWordTokens -= $Count;
          $UnanalyzedWordTypes--;

          # Assign to the word and root a reference to the RP structure.
          $WordsToRPStructures{$ComposedWord} = $RobustnessIndex;
          AddValueToArray( \@{ $RootsToRPStructures{$Root} }, $RobustnessIndex );

        }
      }

      # Update the robustness and frequency of the modified RP structure and
      # link words to roots and patterns.
      $RPStructureRobustness[$RobustnessIndex] = RPStructureRobustness( $RobustnessIndex );
      $RPStructureFrequency[$RobustnessIndex] = RPStructureFrequency( $RobustnessIndex );
      LinkWordsToRootsAndPatterns( $RobustnessIndex );

    }
  }
}




#-------------------------------------------------------------------------------
# FUNCTION: NewRPStructure
#-------------------------------------------------------------------------------
# SYNOPSIS: Constructs a new RP structure, i.e. a hash containing a hash of
# roots and a hash of patterns indexing their token count (within this RP
# structure). Note that RP structures are not really lists of pointers, although
# this is how they are described in the mathematical formalism; the actual
# strings are used as keys in the hashes.
#-------------------------------------------------------------------------------
# ARGUMENTS: References to the lists of roots and the list of patterns of this
# particular RP structure)
#-------------------------------------------------------------------------------
# RETURN VALUE: A reference to the newly constructed RP structure
#-------------------------------------------------------------------------------

sub NewRPStructure {

  my( $RootListRef, $PatternListRef ) = @_;
  my( %RPStructure, $Root, $Pattern, $Count );

  # Loop over roots and patterns and increment their count.

  foreach $Root (@$RootListRef) {
    foreach $Pattern (@$PatternListRef) {

      $Count = $Corpus{ Compose( $Root, $Pattern ) };
      $RPStructure{roots}{$Root} += $Count;
      $RPStructure{patterns}{$Pattern} += $Count;

    }
  }

  return \%RPStructure; # Return a reference to the RP structure.
}



#-------------------------------------------------------------------------------
# FUNCTION: DeleteRPStructure
#-------------------------------------------------------------------------------
# SYNOPSIS: Deletes an RP structure and all references to it.
#-------------------------------------------------------------------------------
# ARGUMENTS: The index of the RP structure to be deleted.
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub DeleteRPStructure {

  my( $RPStructureIndex ) = @_;
  my( $Word, $Root, $Pattern, $Size, $i );

  # Decrement the count of each root and pattern in the RP structure.
  foreach $Root (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} }) {
    foreach $Pattern (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} }) {

      $Word = Compose( $Root, $Pattern );
      $ListOfRoots{$Root} -= $Corpus{$Word};
      $ListOfPatterns{$Pattern} -= $Corpus{$Word};

      # If the corresponding words still points to this RP structure (i.e. if
      # this was not modified somewhere else, e.g. in the collapsing function),
      # mark them as unanalyzed.
      if( $WordsToRPStructures{$Word} == $RPStructureIndex ) {
        $WordsToRPStructures{$Word} = -1;
      }
    }
  }

  undef %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} };
  undef %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} };
  undef %{ $ListOfRPStructures[$RPStructureIndex] };

  # Remove the cell in the @ListOfRPStructures and @RPStructureRobustness arrays.
  RemoveItemFromArray( \@ListOfRPStructures, $RPStructureIndex );

  # Update pointers from words
  foreach $Word (keys %WordsToRPStructures) {
    if( $WordsToRPStructures{$Word} > $RPStructureIndex ) {
      $WordsToRPStructures{$Word}--;
    }
  }

  # Update pointers from roots
  foreach $Root (keys %RootsToRPStructures) {
    RemoveValueFromArray( \@{ $RootsToRPStructures{$Root} }, $RPStructureIndex );

    # If the root points to no more RP structure, delete it.
    if( @{ $RootsToRPStructures{$Root} } == 0 ) {

      delete $ListOfRoots{$Root};
      delete $RootsToRPStructures{$Root};
      delete $RootCreator{$Root};

    }
  }

  # Update pointers from patterns
  foreach $Pattern (keys %PatternsToRPStructures) {
    RemoveValueFromArray( \@{ $PatternsToRPStructures{$Pattern} }, $RPStructureIndex );

    # If the pattern points to no more RP structure, delete it.
    if( @{ $PatternsToRPStructures{$Pattern} } == 0 ) {

      delete $ListOfPatterns{$Pattern};
      delete $PatternsToRPStructures{$Pattern};

    }
  }

}



#-------------------------------------------------------------------------------
# FUNCTION: CollapseRPStructures
#-------------------------------------------------------------------------------
# SYNOPSIS: Moves the whole content of an RP structure to another and updates
# the corresponding pointers.
#-------------------------------------------------------------------------------
# ARGUMENTS: The respective indices of the 'recipient' and 'sender' RP structures
#-------------------------------------------------------------------------------
# RETURN VALUE: none
#-------------------------------------------------------------------------------

sub CollapseRPStructures {

  my( $RPStructureIndex1, $RPStructureIndex2 ) = @_; # 1 = recipient, 2 = sender
  my( $Word, $Root, $Pattern );

  # Increment the count of roots (in the recipient RP structure and in general).
  foreach $Root ( keys %{ ${ $ListOfRPStructures[$RPStructureIndex2] }{roots} } ){

    ${ ${ $ListOfRPStructures[$RPStructureIndex1] }{roots} }{$Root} +=
    ${ ${ $ListOfRPStructures[$RPStructureIndex2] }{roots} }{$Root};

    $ListOfRoots{$Root} +=
    ${ ${ $ListOfRPStructures[$RPStructureIndex2] }{roots} }{$Root};

    # Change pointers from each root to the RP structure.
    AddValueToArray( \@{ $RootsToRPStructures{$Root} }, $RPStructureIndex1 );

    # Change pointers from each word to the RP structure.
    foreach $Pattern
      ( keys %{ ${ $ListOfRPStructures[$RPStructureIndex2] }{patterns} } ) {

      $Word = Compose( $Root, $Pattern );

      $WordsToRPStructures{$Word} = $RPStructureIndex1;

    }
  }

  # Increment the count of patterns (in the RP structure and in general).
  foreach $Pattern
    ( keys %{ ${ $ListOfRPStructures[$RPStructureIndex2] }{patterns} } ){

    ${ ${ $ListOfRPStructures[$RPStructureIndex1] }{patterns} }{$Pattern} +=
    ${ ${ $ListOfRPStructures[$RPStructureIndex2] }{patterns} }{$Pattern};

    $ListOfPatterns{$Pattern} +=
    ${ ${ $ListOfRPStructures[$RPStructureIndex2] }{patterns} }{$Pattern};

    # Change pointers from each pattern to the RP structure.
    AddValueToArray
      ( \@{ $PatternsToRPStructures{$Pattern} }, $RPStructureIndex1 );
  }

}




#-------------------------------------------------------------------------------
# FUNCTION: RPStructureRobustness
#-------------------------------------------------------------------------------
# SYNOPSIS: Calculates the robustness of a given RP structure.
#-------------------------------------------------------------------------------
# ARGUMENTS: The index of the RP structure to be evaluated
#-------------------------------------------------------------------------------
# RETURN VALUE: The robustness of the RP structure
#-------------------------------------------------------------------------------

sub RPStructureRobustness {

  my( $RPStructureIndex ) = @_;
  my( $Root, $Pattern, $TotalLengthOfRoots, $TotalLengthOfPatterns );
  my( $NumOfRoots, $Robustness );

  # Along the lines of Goldsmith 2006, the robustness of an RPStructure is defined
  # as p(n-1) - q, where p and q denote the total length of patterns and
  # roots in the RP structure, and n the number of roots.

  $TotalLengthOfRoots = 0;
  $TotalLengthOfPatterns = 0;

  foreach $Root
    (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} }) {

    $TotalLengthOfRoots += length( $Root );
  }

  foreach $Pattern
    (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} }) {

    $TotalLengthOfPatterns += length( $Pattern );
  }

  $NumOfRoots =
    keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} };

  $Robustness = $TotalLengthOfPatterns * ( $NumOfRoots - 1 )
                - $TotalLengthOfRoots;

  return $Robustness;
}



#-------------------------------------------------------------------------------
# FUNCTION: RPStructureRobustnessOrder
#-------------------------------------------------------------------------------
# SYNOPSIS: Constructs an array where the 1st cell contains the index of the
# most robust RP structure, the 2nd cell contains the index of the second most
# robust RP structure and so on.
#-------------------------------------------------------------------------------
# ARGUMENTS: Reference to the array for storing the ordering.
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub RPStructureRobustnessOrder {

  my( $OrderingRef ) = @_;
  my( @DummyArray, $RPStructureIndex, $NumOfRPStructures );

  undef @$OrderingRef;

  $NumOfRPStructures = int (@RPStructureRobustness);

  # Create a dummy ordering array (1, 2, ...)
  for( $RPStructureIndex = 0;
       $RPStructureIndex < $NumOfRPStructures; $RPStructureIndex++ ) {

       $DummyArray[$RPStructureIndex] = $RPStructureIndex;
  }

  @$OrderingRef =
    sort { $RPStructureRobustness[$b] <=> $RPStructureRobustness[$a] } @DummyArray;

}



#-------------------------------------------------------------------------------
# FUNCTION: RPStructureFrequency
#-------------------------------------------------------------------------------
# SYNOPSIS: Calculates the frequency (types OR tokens) of a given RP structure.
#-------------------------------------------------------------------------------
# ARGUMENTS: The index of the RP structure to be evaluated and optionally the
# string 'types' indicating that the counting mode is not the default (tokens).
#-------------------------------------------------------------------------------
# RETURN VALUE: The frequency of the RP structure
#-------------------------------------------------------------------------------

sub RPStructureFrequency {

  my( $RPStructureIndex, $Mode ) = @_;
  my( $Pattern, $Root, $Frequency );

  $Frequency = 0;

  unless( defined $Mode ) { $Mode = 'tokens'; }

  foreach $Pattern
    (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} }) {

    if( $Mode ne 'types' ) {
      $Frequency += ${ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} }{$Pattern};
    }
    else {
      $Frequency += int (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} });
    }

  }

  return $Frequency;
}



#-------------------------------------------------------------------------------
# FUNCTION: LinkWordsToRootsAndPatterns
#-------------------------------------------------------------------------------
# SYNOPSIS: Create pointers to the roots and patterns of a given RP structure
# from the corresponding words.
#-------------------------------------------------------------------------------
# ARGUMENTS: The index of the structure whose roots and patterns must be linked
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub LinkWordsToRootsAndPatterns {

  my( $RPStructureIndex ) = @_;
  my( $Word, $Root, $Pattern );

  foreach $Root
    (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} }) {

    foreach $Pattern
      (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} }) {

      $Word = Compose( $Root, $Pattern );

      $WordsToRoots{ $Word } = $Root;
      $WordsToPatterns{ $Word } = $Pattern;
    }
  }

}



#-------------------------------------------------------------------------------
# FUNCTION: Decompose
#-------------------------------------------------------------------------------
# SYNOPSIS: Analyzes a word into a root and a pattern.
#-------------------------------------------------------------------------------
# ARGUMENTS: The word and a reference to a hash for storing the analysis
#-------------------------------------------------------------------------------
# RETURN VALUE: None.
#-------------------------------------------------------------------------------

sub Decompose {

  my( $Word, $AnalysisRef ) = @_;
  my( $Length, $Phoneme, $Vowel, $RootCategory, $i );

  $Length = length( $Word );

  # Check if patterns are associated with consonants or vowels
  # in the parameters.

  if( $Parameters{PhonemeCategoryInPatterns} eq 'vowels' ) {
    $RootCategory = \%Consonants;
  }
  else { $RootCategory = \%Vowels; }

  $$AnalysisRef{root} = "";
  $$AnalysisRef{pattern} = "";

  # For each phoneme
  for( $i = 0; $i < $Length; $i++ ) {

    $Phoneme = substr( $Word, $i, 1 );

    # If it is the category selected for roots, add it to the root and
    # add a C slot to the pattern.
    if( defined $$RootCategory{$Phoneme} ) {
      $$AnalysisRef{root} .= $Phoneme;
      $$AnalysisRef{pattern} .= $Parameters{InsertionSlotSymbol};
    }

    # Else add it to the pattern.
    else { $$AnalysisRef{pattern} .= $Phoneme; }
  }

  return;
}



#-------------------------------------------------------------------------------
# FUNCTION: Compose
#-------------------------------------------------------------------------------
# SYNOPSIS: Turns a pair of root and pattern into an actual word.
#-------------------------------------------------------------------------------
# ARGUMENTS: The root and pattern
#-------------------------------------------------------------------------------
# RETURN VALUE: The word
#-------------------------------------------------------------------------------

sub Compose {

  my( $Root, $Pattern ) = @_;
  my( $Word, $Length, $Category, $i, $RootCategoryIndex );

  $RootCategoryIndex = 0;

  $Length = length( $Pattern );

  # For each symbol in the pattern
  for( $i = 0; $i < $Length; $i++ ) {

    $Category = substr( $Pattern, $i, 1 );

    # If it is an insertion slot, look it up in the root and add it
    # to the word.
    if( $Category eq $Parameters{InsertionSlotSymbol} ) {
      $Word .= substr( $Root, $RootCategoryIndex++, 1 );
    }

    # Else add it to the word.
    else { $Word .= $Category; }
  }

  return( $Word );

}



#-------------------------------------------------------------------------------
# FUNCTION: Sukhotin
#-------------------------------------------------------------------------------
# SYNOPSIS: Classifies phonemes into consonants and vowels, using Sukhotin's
# algorithm (Sukhotin, 1962). Consonants are stored into the (global)
# %Consonants hash, and vowels into %Vowels.
#-------------------------------------------------------------------------------
# ARGUMENTS: none
#-------------------------------------------------------------------------------
# RETURN VALUE: none
#-------------------------------------------------------------------------------

sub Sukhotin {

  my( $Word, $Length, $TableRow, $LeftPhoneme, $RightPhoneme, %UnigramCount );
  my( %BigramCount, $Phoneme, $NumOfPhonemes, $MaxPhoneme, $VowelsFound );
  my( $WordCount, $i );

  # First count unigrams and (unordered) bigrams frequencies.

  foreach $Word ( keys %Corpus ) {

    $WordCount = $Corpus{$Word};

    $Length = length( $Word );

    for( $i = 0; $i < $Length-1; $i++ ) { # Ignore repeated symbols.

      $LeftPhoneme = substr( $Word, $i, 1 );
      $RightPhoneme = substr( $Word, $i+1, 1 );

      if( $RightPhoneme ne $LeftPhoneme ) {

        $BigramCount{$LeftPhoneme}{$RightPhoneme} += $WordCount;
        $BigramCount{$RightPhoneme}{$LeftPhoneme} += $WordCount;

        $UnigramCount{$LeftPhoneme} += $WordCount;
        $UnigramCount{$RightPhoneme} += $WordCount;
      }
    }
  }

  $NumOfPhonemes = int (keys %UnigramCount);

  $VowelsFound = 0;

  # For each phoneme (i.e. row in the table of unordered pairs of phonemes)
  for( $TableRow = 0; $TableRow < $NumOfPhonemes; $TableRow++ ) {

    # Find the phoneme with maximal frequency (in bigrams).
    $MaxPhoneme = MaxSum( \%UnigramCount, \$VowelsFound );

    # If all vowels have been found, exit the loop.
    if( $VowelsFound ) { last; }

    # Else label the phoneme with maximal sum as a vowel.
    $Vowels{$MaxPhoneme} = 1;

    delete $UnigramCount{$MaxPhoneme}; # Remove it from the table.

    # Update the sum for each phoneme (i.e. column in the table).
    foreach $Phoneme (keys %UnigramCount) {
      if( exists $BigramCount{$MaxPhoneme}{$Phoneme} ) {
        $UnigramCount{$Phoneme} -= 2 * $BigramCount{$MaxPhoneme}{$Phoneme};
      }
    }

  }

  # Label remaining phonemes as consonants.
  foreach $Phoneme (keys %UnigramCount) {
    $Consonants{$Phoneme} = 1;
  }

}



#-------------------------------------------------------------------------------
# FUNCTION: MaxSum
#-------------------------------------------------------------------------------
# SYNOPSIS: During Sukhotin's algorithm, finds the phoneme (i.e. table row)
# with the largest sum, and indicates when all vowels have been found.
#-------------------------------------------------------------------------------
# ARGUMENTS: References to the hash containing the sum for each phoneme and to
#            the flag signalling that all vowels were found.
#-------------------------------------------------------------------------------
# RETURN VALUE: The phoneme with maximal sum.
#-------------------------------------------------------------------------------

sub MaxSum {

  my( $UnigramCountRef, $VowelsFoundRef ) = @_;
  my( $Phoneme, $Max, $MaxPhoneme );

  $Max = -1;

  foreach $Phoneme (keys %$UnigramCountRef) {

    if( $$UnigramCountRef{$Phoneme} > $Max ) {

      $Max = $$UnigramCountRef{$Phoneme};
      $MaxPhoneme = $Phoneme;

    }
  }

  # If no positive sum remains, indicate that all vowels have been found.
  if( $Max <= 0 ) { $$VowelsFoundRef = 1; }

  return $MaxPhoneme;
}



#-------------------------------------------------------------------------------
# FUNCTION: AddValueToArray
#-------------------------------------------------------------------------------
# SYNOPSIS: Adds a specific value to an array if it is not already in it.
#-------------------------------------------------------------------------------
# ARGUMENTS: A reference to the array and the value to be added
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub AddValueToArray {

  my( $ArrayRef, $NewValue ) = @_;
  my( $Found, $Value );

  $Found = 0;

  foreach $Value (@$ArrayRef) {
    if( $Value == $NewValue ) { $Found = 1; last; }
  }

  unless( $Found ) { push @$ArrayRef, $NewValue }

}



#-------------------------------------------------------------------------------
# FUNCTION: RemoveValueFromArray
#-------------------------------------------------------------------------------
# SYNOPSIS: Removes a specific value from an array and decreases all larger
# values by 1.
#-------------------------------------------------------------------------------
# ARGUMENTS: A reference to the array and the value to be removed
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub RemoveValueFromArray {

  my( $ArrayRef, $ValueToRemove ) = @_;
  my( $Size, $Found, $ItemToRemove, $i );

  $Found = 0;

  $Size = int (@$ArrayRef);

  for( $i = 0; $i < $Size; $i++ ) {

    if( $$ArrayRef[$i] == $ValueToRemove ) { $Found = 1; $ItemToRemove = $i }
    elsif( $$ArrayRef[$i] > $ValueToRemove ) { $$ArrayRef[$i]--; }

  }

  if( $Found ) { RemoveItemFromArray( $ArrayRef, $ItemToRemove ); }

}



#-------------------------------------------------------------------------------
# FUNCTION: RemoveItemFromArray
#-------------------------------------------------------------------------------
# SYNOPSIS: Removes an item from an array.
#-------------------------------------------------------------------------------
# ARGUMENTS: A reference to the array and the index of the item to be removed
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub RemoveItemFromArray {

  my( $ArrayRef, $ItemToRemove ) = @_;
  my( $Item, $Size, $Found, $i );

  $Found = 0;

  $Size = int (@$ArrayRef);

  for( $i = 0; $i < $Size; $i++ ) {

    if( $i > $ItemToRemove ) {
      $$ArrayRef[$i-1] = $$ArrayRef[$i];
    }

  }

  pop @$ArrayRef;

}



#-------------------------------------------------------------------------------
# FUNCTION: CompareArrays
#-------------------------------------------------------------------------------
# SYNOPSIS: Checks if 2 (unordered) arrays contain exactly the same elements.
#-------------------------------------------------------------------------------
# ARGUMENTS: A reference to each array
#-------------------------------------------------------------------------------
# RETURN VALUE: 1 if the arrays are identical, 0 otherwise
#-------------------------------------------------------------------------------

sub CompareArrays {

  my( $ArrayRef1, $ArrayRef2 ) = @_;
  my( @SortedArray1, @SortedArray2, $Result, $Size1, $Size2, $i );

  $Result = 1;

  $Size1 = int (@$ArrayRef1);
  $Size2 = int (@$ArrayRef2);

  if( $Size1 != $Size2 ) { return 0; }

  @SortedArray1 = sort @$ArrayRef1;
  @SortedArray2 = sort @$ArrayRef2;

  for( $i = 0; $i < $Size1; $i++ ) { # NB: we assume that STRINGS are compared.
    if( $SortedArray1[$i] ne $SortedArray2[$i] ) { $Result = 0; last; }
  }

  return $Result;

}





#-------------------------------------------------------------------------------
# FUNCTION: DeltaDescriptionLength
#-------------------------------------------------------------------------------
# SYNOPSIS: Computes the difference in description length resulting from the
# addition of a set of new roots to an existing RP structure.
#-------------------------------------------------------------------------------
# ARGUMENTS: A reference to the list of new roots.
#-------------------------------------------------------------------------------
# RETURN VALUE: The difference in description length (DL)
#-------------------------------------------------------------------------------

sub DeltaDescriptionLength {

  my( $RPStructureIndex, $NewRootListRef ) = @_;
  my( $Word, $Root, $Pattern, $ComposedWord, $DeltaDL, $Sum1, $Sum2, $Sum3 );
  my( $NumOfNewRootTypes, $NumOfNewWordTypes, $NumOfRootTypesInRPStructure );
  my( $RootListRef, $PatternListRef, %PatternCount, @TempArray );
  my( $NumOfNewWordTokens );

  $PatternListRef = \%{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} };
  $RootListRef = \%{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} };

  $DeltaDL = 0;

  $NumOfNewRootTypes = int keys %$NewRootListRef;

  $NumOfNewWordTypes = $NumOfNewRootTypes * int keys %$PatternListRef;

  $NumOfRootTypesInRPStructure = int keys %$RootListRef;

  #------ Changes in model DL

  # Change in structural cost of the general list of roots:
  $DeltaDL += ListCost( ( ( int( keys %ListOfRoots ) + $UnanalyzedWordTypes )
                          - ( $NumOfNewWordTypes - $NumOfNewRootTypes ) )
                        / ( int( keys %ListOfRoots ) + $UnanalyzedWordTypes ) ); #/

  # Change in cost of the list of pointers to roots in the RP structure:
  $DeltaDL += ListCost( ( $NumOfRootTypesInRPStructure + $NumOfNewRootTypes )
                        / $NumOfRootTypesInRPStructure );                     #/

  # Change in cost of the list of pointers to roots in the NULL RP structure:
  $DeltaDL += ListCost( ( $UnanalyzedWordTypes - $NumOfNewWordTypes )
                        / $UnanalyzedWordTypes );                             #/

  # Addition and suppression of roots ($Sum1), addition of pointers to roots in
  # RP structure ($Sum2), and suppression of pointers to roots in NULL RP
  # structure ($Sum3) :

  @TempArray = keys %$NewRootListRef;
  $Sum1 = $NumOfNewRootTypes * length $TempArray[0];

  $Sum2 = $NumOfNewRootTypes * Log2( $NumWords );

  $Sum3 = -$NumOfNewWordTypes * Log2( $NumWords );

  foreach $Root (keys %$NewRootListRef) {

    $Sum2 -= Log2( $$NewRootListRef{$Root} );

    foreach $Pattern (keys %$PatternListRef) {

      $ComposedWord = Compose( $Root, $Pattern );

      $Sum1 -= length $ComposedWord;

      $Sum3 += Log2( $Corpus{$ComposedWord} );

      $PatternCount{$Pattern} += $Corpus{$ComposedWord};

      $NumOfNewWordTokens += $Corpus{$ComposedWord};
    }

  }

  $DeltaDL += Log2( int (keys %Consonants) + int (keys %Vowels) ) * $Sum1;
  $DeltaDL += $Sum2;
  $DeltaDL += $Sum3;

  # Change in cost of pointers to patterns in RP structures

  $Sum1 = 0;

  foreach $Pattern (keys %$PatternListRef) {

    $Sum1 -= int( @{ $PatternsToRPStructures{$Pattern} })
           * Log2( 1 + ( $PatternCount{$Pattern} / $ListOfPatterns{ $Pattern } ) );

  }

  $DeltaDL += $Sum1;

  # Change in cost of pointer to NULL pattern
  $DeltaDL -= Log2( 1 - ( $NumOfNewWordTokens / $UnanalyzedWordTokens ) );

  #------ Changes in compressed corpus DL

  # Change in cost of new analyzed words ($Sum1) and words in RP structure ($Sum2)
  $Sum1 = 0;
  $Sum2 = 0;

  foreach $Pattern (keys %$PatternListRef) {

    foreach $Root (keys %$NewRootListRef) {

      $ComposedWord = Compose( $Root, $Pattern );

      $Sum1 += $Corpus{$ComposedWord} * Log2 ( ( $Corpus{$ComposedWord}
                       * ( $RPStructureFrequency[$RPStructureIndex] + $NumOfNewWordTokens ) )
                     / ( $$NewRootListRef{$Root}                              #/
                       * ( $PatternCount{$Pattern} + $$PatternListRef{$Pattern} ) ) );

    }

    foreach $Root (keys %$RootListRef) {

      $ComposedWord = Compose( $Root, $Pattern );

      $Sum2 += $Corpus{$ComposedWord} * Log2 ( ( ( $RPStructureFrequency[$RPStructureIndex] + $NumOfNewWordTokens )
                        * $$PatternListRef{$Pattern} )
                      / ( $RPStructureFrequency[$RPStructureIndex]            #/
                        * ( $PatternCount{$Pattern} + $$PatternListRef{$Pattern} ) ) );
    }

  }

  $DeltaDL += $Sum1;
  $DeltaDL += $Sum2;

  return $DeltaDL;
}





#-------------------------------------------------------------------------------
# FUNCTION: RPDescriptionLength
#-------------------------------------------------------------------------------
# SYNOPSIS: Computes the description length of the RP analysis. Stores the
# results in the %RPDescriptionLength hash.
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: None
#-------------------------------------------------------------------------------

sub RPDescriptionLength {

  my( $Word, $Root, $Pattern, $PatternCategory, $RPStructureIndex, $DLWord );
  my( $NumPhon, $NumPhonPatt, $DLRoots, $DLPatterns, $DLRPStructures, $DLCorpus );
  my( $NumRootsInRPStructure, $NumPatternsInRPStructure, $NumOfRPStructures );
  my( $DLNullRPAnalysis );

  if( $Parameters{PhonemeCategoryInPatterns} eq 'vowels' ) {
    $PatternCategory = \%Vowels;
  }
  else { $PatternCategory = \%Consonants; }

  # Num. of distinct phonemes, and distinct symbols in patterns.
  $NumPhon = int (keys %Consonants) + int (keys %Vowels);
  $NumPhonPatt = int (keys %$PatternCategory) + 1;

  #--------- List of roots

  $DLRoots = 0;

  foreach $Root (keys %ListOfRoots) {
    $DLRoots += length( $Root );
  }

  foreach $Word (keys %Corpus) {

    if( defined $WordsToRoots{$Word} == 0 ) { # Add unanalyzed words
      $DLRoots += length( $Word );
    }

  }

  $DLRoots *= Log2( $NumPhon );

  $DLRoots += ListCost( int ( keys %ListOfRoots ) + $UnanalyzedWordTypes );

  #--------- List of patterns

  $DLPatterns = 0;

  foreach $Pattern (keys %ListOfPatterns) {
    $DLPatterns += length( $Pattern );
  }

  $DLPatterns *= Log2( $NumPhonPatt );

  $DLPatterns += ListCost( int ( keys %ListOfPatterns ) + 1 );

  #--------- List of RP structures

  $NumOfRPStructures = int ( @ListOfRPStructures );

  if( $UnanalyzedWordTokens > 0 ) {
    $DLRPStructures = ListCost( $NumOfRPStructures + 1 );
  }
  else{ $DLRPStructures = ListCost( $NumOfRPStructures ); }
  for( $RPStructureIndex = 0; $RPStructureIndex < $NumOfRPStructures; $RPStructureIndex++ ) {

    $NumRootsInRPStructure =
        int (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} } );
    $NumPatternsInRPStructure =
        int (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} } );

    $DLRPStructures += ListCost( $NumRootsInRPStructure );
    $DLRPStructures += ListCost( $NumPatternsInRPStructure );

    $DLRPStructures += ( $NumRootsInRPStructure + $NumPatternsInRPStructure ) * Log2( $NumWords );

    foreach $Root (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{roots} }) {
      $DLRPStructures -= Log2( $ListOfRoots{$Root} );
    }

    foreach $Pattern (keys %{ ${ $ListOfRPStructures[$RPStructureIndex] }{patterns} }) {
      $DLRPStructures -= Log2( $ListOfPatterns{$Pattern} );
    }

  }

  if( $UnanalyzedWordTokens > 0 ) { # Add Null RP structure

    $DLRPStructures += ListCost( $UnanalyzedWordTypes );

    $DLRPStructures += ( $UnanalyzedWordTypes + 1 ) * Log2( $NumWords );

    foreach $Word (keys %Corpus) {

      if( defined $WordsToRoots{$Word} == 0 ) {
        $DLRPStructures -= Log2( $Corpus{$Word} );
      }

    }

    $DLRPStructures -= Log2( $UnanalyzedWordTokens );

  }

  #--------- Compressed corpus

  $DLCorpus = $NumWords * Log2( $NumWords );

  foreach $Word (keys %Corpus) {

    if( defined $WordsToRoots{$Word} ) { # Analyzed words

      $DLWord =
        Log2( $RPStructureFrequency[ $WordsToRPStructures{ $Word } ]
              / ( $ListOfRoots{ $WordsToRoots{ $Word } }                #/
                  * ${ ${ $ListOfRPStructures[ $WordsToRPStructures{ $Word } ] }
                        {patterns} }{$WordsToPatterns{ $Word }} ) );

      $DLCorpus += $Corpus{$Word} * $DLWord;

    }
    else { # Unanalyzed words

      $DLWord = Log2( $Corpus{ $Word } );

      $DLCorpus -= $Corpus{$Word} * $DLWord;

    }

  }

  #--------- Null RP Analysis

  $DLNullRPAnalysis = 2 * ListCost( int (keys %Corpus) );

  $DLNullRPAnalysis += ( $NumWords + int (keys %Corpus) ) * Log2( $NumWords );

  foreach $Word (keys %Corpus) {

    $DLNullRPAnalysis += ( length( $Word ) * Log2( $NumPhon ) )
                       - ( $Corpus{$Word}+1 ) * Log2( $Corpus{$Word} );

  }

  $RPDescriptionLength{roots} = $DLRoots;
  $RPDescriptionLength{patterns} = $DLPatterns;
  $RPDescriptionLength{rpstructures} = $DLRPStructures;
  $RPDescriptionLength{corpus} = $DLCorpus;
  $RPDescriptionLength{total} = $DLRoots + $DLPatterns + $DLRPStructures + $DLCorpus;
  $RPDescriptionLength{nullrpanalysis} = $DLNullRPAnalysis;
  $RPDescriptionLength{deltarpanalysis} = $DLNullRPAnalysis - $RPDescriptionLength{total};
  $RPDescriptionLength{relativedeltarpanalysis} =
              100 * ($RPDescriptionLength{deltarpanalysis} / $DLNullRPAnalysis);

}




#-------------------------------------------------------------------------------
# FUNCTION: ListCost
#-------------------------------------------------------------------------------
# SYNOPSIS: Calculates the cost in bits of a list structure of given length.
#-------------------------------------------------------------------------------
# ARGUMENTS: The length of the list
#-------------------------------------------------------------------------------
# RETURN VALUE: The cost of the list structure
#-------------------------------------------------------------------------------

sub ListCost {

  my( $Length ) = @_;

  if( $Length > 0 ) { return Log2( $Length ); }
  else { return 0; }

}



#-------------------------------------------------------------------------------
# FUNCTION: Log2
#-------------------------------------------------------------------------------
# SYNOPSIS: Returns the base 2 log of the argument (or 0).
#-------------------------------------------------------------------------------
# ARGUMENTS: A positive number
#-------------------------------------------------------------------------------
# RETURN VALUE: Its base 2 log
#-------------------------------------------------------------------------------

sub Log2 {

  my( $Number ) = @_;

  if( $Number != 0 ) { return log( $Number ) / log( 2 ); }
  else{ return 0; }

}



#-------------------------------------------------------------------------------
# FUNCTION: NumberOfPhonemesInPattern
#-------------------------------------------------------------------------------
# SYNOPSIS: Calculates the number of vowels in a pattern.
#-------------------------------------------------------------------------------
# ARGUMENTS: The pattern
#-------------------------------------------------------------------------------
# RETURN VALUE: The number of vowel
#-------------------------------------------------------------------------------

sub NumberOfPhonemesInPattern {

  my( $Pattern ) = @_;

  $Pattern =~ s/$Parameters{InsertionSlotSymbol}//g;

  return length( $Pattern );

}




#-------------------------------------------------------------------------------
# FUNCTION: OutputResults
#-------------------------------------------------------------------------------
# SYNOPSIS: Outputs a summary of the induced morphology to a file (or STDOUT).
#-------------------------------------------------------------------------------
# ARGUMENTS: None
#-------------------------------------------------------------------------------
# RETURN VALUE: 1 if the file was written correctly, 0 otherwise
#-------------------------------------------------------------------------------

sub OutputResults {

  my( $Word, $String, @RPStructureRobustnessOrder, $Index, $Root, $Pattern, $Item );
  my( $RPStructures, $Output, $RPStructureIndex, $HLine, $Phoneme );
  my( $DL, $DL2, $NumOfRPStructures, $NumPhon, $NumPhonPatt, $PatternCategory );
  my( $NumRootsInRPStructure, $NumPatternsInRPStructure );

  # Determine if the output should be sent to a file or to the console (STDOUT).

  if( defined $Parameters{ OutputFileName } ) {
    open( $Output, ">$Parameters{ OutputFileName }" ) || return 0;
  }
  else {
    open( $Output, ">&STDOUT" );
  }

  # Determine whether patterns are associated with consonants or vowels
  # in the parameters.

  if( $Parameters{ PhonemeCategoryInPatterns } eq 'vowels' ) {
    $PatternCategory = \%Vowels;
  }
  else { $PatternCategory = \%Consonants; }

  # Num. of distinct phonemes, and distinct symbols in patterns.

  $NumPhon = int (keys %Consonants) + int (keys %Vowels);
  $NumPhonPatt = int (keys %$PatternCategory) + 1;

  $HLine = "-------------------------------------------------------------\n";

  # Main header

  print $Output "\n$HLine* ROOT-AND-PATTERN ANALYSIS FOR FILE ";
  print $Output $Parameters{ InputFileName }."\n$HLine";

  # Description length

  print $Output "\n$HLine"."* DESCRIPTION LENGTH ANALYSIS\n".$HLine;

  print $Output "List of roots:\t";
  print $Output sprintf( "%.2f", $RPDescriptionLength{ roots } )." bits\n";

  print $Output "List of patterns:\t";
  print $Output sprintf( "%.2f", $RPDescriptionLength{ patterns } )." bits\n";

  print $Output "List of RP structures:\t";
  print $Output
    sprintf( "%.2f", $RPDescriptionLength{ rpstructures } )." bits\n";

  print $Output "Compressed corpus:\t";
  print $Output sprintf( "%.2f", $RPDescriptionLength{ corpus } )." bits\n";

  print $Output "Total:\t";
  print $Output sprintf( "%.2f", $RPDescriptionLength{ total } )." bits\n";

  print $Output "Delta:\t";
  print $Output
    sprintf( "%.2f", $RPDescriptionLength{ deltarpanalysis } )." bits\n";

  print $Output "Relative delta:\t";
  print $Output
    sprintf( "%.2f", $RPDescriptionLength{ relativedeltarpanalysis } )."%\n";

  # Words

  print $Output "\n$HLine"."* WORDS: ";
  print $Output int(keys %Corpus)." types and $NumWords tokens\n";

  print $Output "* DL: ";
  print $Output sprintf("%.2f", $RPDescriptionLength{ corpus })." bits\n";
  print $Output $HLine;

  print $Output "WORD\tCOUNT\tRP STRUCTURE\tDL\n";

  foreach $Word (sort keys %Corpus) {

    if( defined $WordsToRoots{ $Word } ) {

      # Analyzed words

      $DL =
        Log2( $NumWords * $RPStructureFrequency[ $WordsToRPStructures{ $Word } ]
              / ( $ListOfRoots{ $WordsToRoots{ $Word } }
                  * ${ ${ $ListOfRPStructures[ $WordsToRPStructures{ $Word } ] }
                        {patterns} }{$WordsToPatterns{ $Word }} ) );
    }
    else {

      # Unanalyzed words

      $DL = Log2( $NumWords / $Corpus{ $Word } );

    }
    print $Output "$Word\t$Corpus{ $Word }\t$WordsToRPStructures{ $Word }\t";
    print $Output sprintf( "%.2f", $DL )."\n";

  }

  # Phonological analysis

  print $Output "\n$HLine"."* CONSONANTS AND VOWELS\n".$HLine;

  $String = "C = {";
  foreach $Phoneme (sort keys %Consonants) {
    $String .= "$Phoneme,";
  }
  chop $String;

  $String .= "}\nV = {";
  foreach $Phoneme (sort keys %Vowels) {
    $String .= "$Phoneme,";
  }
  chop $String;

  print $Output "$String}\n";

  # Roots

  print $Output "\n$HLine"."* ROOTS: ";
  print $Output (int(keys %ListOfRoots)+$UnanalyzedWordTypes);
  print $Output " types and $NumWords tokens\n";

  print $Output "* DL: ";
  print $Output sprintf( "%.2f", $RPDescriptionLength{ roots }
                                 - ListCost( int ( keys %ListOfRoots )
                                             + $UnanalyzedWordTypes ) );
  print $Output " bits";

  print $Output " (+ ";
  print $Output sprintf( "%.2f", ListCost( int ( keys %ListOfRoots )
                                           + $UnanalyzedWordTypes ) ).")\n";
  print $Output $HLine;

  print $Output "ROOT\tCOUNT\tRP STRUCTURE\tDL";

  if( $UnanalyzedWordTokens < $NumWords ) { print $Output "\tCreator"; }

  print $Output "\n";

  foreach $Root (sort keys %ListOfRoots) {

    print $Output "$Root\t$ListOfRoots{ $Root }\t";

# CODE FOR DISPLAYING SEVERAL RP STRUCTURES PER ROOT
#    $String = "\t{";
#    foreach $Item (sort { $a <=> $b } @{ $RootsToRPStructures{$Root} }) {
#      $String .= "$Item,";
#    }
#    chop $String;
#    print $Output "$String}\n";

    print $Output "${ $RootsToRPStructures{ $Root } }[0]\t";
    print $Output sprintf( "%.2f", length( $Root )*Log2( $NumPhon ) );
    if( defined $RootCreator{ $Root } ) {
      print $Output "\t$RootCreator{ $Root }";
    }

    print $Output "\n";

  }

  # Add unanalyzed words

  foreach $Word (sort keys %Corpus) {

    if( defined $WordsToRoots{ $Word } == 0 ) {
      print $Output "$Word\t$Corpus{ $Word }\t-1\t";
      print $Output sprintf( "%.2f",length( $Word )*Log2( $NumPhon ) )."\n";
    }

  }

  # Patterns

  print $Output "\n$HLine"."* PATTERNS: ";
  print $Output int (keys %ListOfPatterns)." types and $NumWords tokens\n";

  print $Output "* DL: ";

  if( $UnanalyzedWordTokens > 0 ) {

    print $Output
      sprintf( "%.2f", $RPDescriptionLength{ patterns }
                       - ListCost( int (keys %ListOfPatterns) + 1 ) );
    print $Output " bits (+ ";
    print $Output
      sprintf( "%.2f", ListCost( int ( keys %ListOfPatterns ) + 1 ) ).")\n";

  }
  else {

    print $Output
      sprintf( "%.2f", $RPDescriptionLength{ patterns }
                       - ListCost( int ( keys %ListOfPatterns ) ) );
    print $Output " bits (+ ";
    print $Output
      sprintf( "%.2f", ListCost( int ( keys %ListOfPatterns  ) ) ).")\n";

  }

  print $Output $HLine;
  print $Output "PATTERN\tCOUNT\tRP STRUCTURES\tDL\n";

  if( $UnanalyzedWordTokens > 0  ) {
    print $Output "NULL\t$UnanalyzedWordTokens\t{-1}\t0\n";
  }

  foreach $Pattern (sort keys %ListOfPatterns) {

    print $Output "$Pattern\t$ListOfPatterns{ $Pattern }";

    $String = "\t{";

    foreach $Item
      (sort { $a <=> $b } @{ $PatternsToRPStructures{ $Pattern } })
    {
      $String .= "$Item,";
    }
    chop $String;

    print $Output "$String}\t";
    print $Output
      sprintf( "%.2f",length( $Pattern )*Log2( $NumPhonPatt ) )."\n";

  }

  # RP structures

  print $Output "\n$HLine"."* RP STRUCTURES: ";
  print $Output int (@ListOfRPStructures+1)." types and $NumWords tokens\n";

  print $Output "* DL: ";
  print $Output
    sprintf( "%.2f", $RPDescriptionLength{ rpstructures }
                     - ListCost( int (@ListOfRPStructures+1) ) );
  print $Output " bits (+ ";
  print $Output
    sprintf( "%.2f", ListCost( int (@ListOfRPStructures+1) ) ).")\n";

  $NumOfRPStructures = int (@ListOfRPStructures);

  RPStructureRobustnessOrder( \@RPStructureRobustnessOrder );

  for( $RPStructureIndex = 0; $RPStructureIndex < $NumOfRPStructures;
    $RPStructureIndex++ )
  {
    print $Output $HLine."INDEX\tCOUNT\tROBUSTNESS\tDL ROOTS\tDL PATTERNS\n";

    $Index = $RPStructureRobustnessOrder[ $RPStructureIndex ];

    $String = "\t-------------\n\tROOT\tCOUNT\tDL\n";

    $DL = 0;
    $DL2 = 0;

    $NumRootsInRPStructure =
      int (keys %{ ${ $ListOfRPStructures[ $Index ] }{ roots } } );

    $NumPatternsInRPStructure =
      int (keys %{ ${ $ListOfRPStructures[ $Index ] }{ patterns } } );

    foreach $Root
      (sort keys %{ ${ $ListOfRPStructures[ $Index ] }{ roots } })
    {
      $String .= "\t$Root\t";
      $String .= ${ ${ $ListOfRPStructures[ $Index ] }{ roots } }{ $Root };
      $String .= "\t";
      $String .=
        sprintf( "%.2f", Log2( $NumWords / $ListOfRoots{ $Root } ) )."\n";

      $DL += Log2( $NumWords / $ListOfRoots{ $Root } );
    }

    $String .= "\t-------------\n\tPATTERN\tCOUNT\n";

    foreach $Pattern
      (sort keys %{ ${ $ListOfRPStructures[ $Index ] }{ patterns } })
    {
      $String .= "\t$Pattern\t";
      $String .=
        ${ ${ $ListOfRPStructures[ $Index ] }{ patterns } }{ $Pattern };
      $String .= "\t";
      $String .=
        sprintf( "%.2f", Log2( $NumWords / $ListOfPatterns{ $Pattern } ) )."\n";

      $DL2 += Log2( $NumWords / $ListOfPatterns{ $Pattern } );
    }

    chop $String;

    print $Output "$Index\t$RPStructureFrequency[ $Index ]\t";
    print $Output "$RPStructureRobustness[ $Index ]\t";
    print $Output sprintf( "%.2f", $DL )." bits (+";
    print $Output sprintf( "%.2f", ListCost( $NumRootsInRPStructure ) ).")\t";
    print $Output sprintf( "%.2f", $DL2 )." bits (+";
    print $Output
      sprintf( "%.2f", ListCost( $NumPatternsInRPStructure ) ).")\n";

    print $Output "$String\n";

  }

  # NULL RP structure

  print $Output $HLine."INDEX\tCOUNT\tROBUSTNESS\tDL ROOTS\tDL PATTERNS\n";

  $String = "\t-------------\n\tROOT\tCOUNT\n";

  $DL = 0;
  $DL2 = 0;

  foreach $Word (sort keys %Corpus) {

    if( defined $WordsToRoots{ $Word } == 0 ) {

      $String .= "\t$Word\t$Corpus{ $Word }\t";
      $String .= sprintf( "%.2f", Log2( $NumWords / $Corpus{ $Word } ) )."\n";

      $DL += Log2( $NumWords / $Corpus{ $Word } );

    }

  }

  $String .= "\t-------------\n\tPATTERN\tCOUNT\n";
  $String .= "\tNULL\t$UnanalyzedWordTokens\tDL\t";

  if( $UnanalyzedWordTokens > 0 ) {
    $String .=
      sprintf( "%.2f", Log2( $NumWords / $UnanalyzedWordTokens ) )."\n";

    $DL2 += Log2( $NumWords / $UnanalyzedWordTokens );
  }

  else {
    $String .= "0\n";
  }

  print $Output "-1\t$UnanalyzedWordTokens\t0\t";
  print $Output sprintf( "%.2f", $DL )." bits (+";
  print $Output sprintf( "%.2f", ListCost( $UnanalyzedWordTypes ) ).")\t";
  print $Output sprintf( "%.2f", $DL2 )." bits\n";

  print $Output "$String";

  close $Output;

  return 1;

}



