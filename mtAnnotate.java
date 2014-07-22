
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Map.Entry;
import java.io.*;

import org.biojava3.core.sequence.*;
import org.biojava3.core.sequence.io.*;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.features.TextFeature;
import org.biojava3.core.sequence.loader.GenbankProxySequenceReader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.util.StringManipulationHelper;
import org.biojava3.alignment.*;
import org.biojava3.alignment.template.*;

public class mtAnnotate {

   public static SubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix
         <NucleotideCompound>(AmbiguityDNACompoundSet.getDNACompoundSet(),
               (short)1,(short)-3);
   public static SimpleGapPenalty gapP = new SimpleGapPenalty((short)5,(short)2);

   public static String getDate() {
      Date sysdate = Calendar.getInstance().getTime();

      String default_date = new SimpleDateFormat("dd-MMM-yyyy")
      .format(sysdate);
      return default_date;
   }
   public static void printHelp(){
      System.out.println(
            "Usage: mtAnnotate [options] GenBankFile FASTAFile\n"
                  + "\n"
                  + "This program takes two closely related DNA sequences, and attempts to align\n"
                  + "them using the Needleman-Wunsch global alignment algorithm. The alignment is\n"
                  + "then used to infer annotations in the inputted FASTA sequence. Basic sequence\n"
                  + "data structures and manipulations, along with the alignment algorithm are\n"
                  + "from the BioJava3 library, with some modifications.\n"
                  + "\n"
                  +"Inputs:  GenBank file containing annotated DNA sequence\n"
                  + "        FASTA file containing closely related raw DNA sequence\n"
                  + "Outputs: GenBank file containing DNA sequence from inputted FASTA file\n"
                  + "         with annotations inferred from alignment.\n"
                  + "        FASTA file containing DNA sequences of annotated features in\n"
                  + "         the inputted FASTA sequence.\n"
                  + "\n"
                  + "Options:\n"
                  + "      -t feature_type...\n"
                  + "         Specifies GenBank feature types for which FASTA sequence\n"
                  + "         outputs are desired, if not specified, FASTA sequences will\n"
                  + "         will be outputted for all features.\n"
                  + "      -g filename\n"
                  + "         Specifies the name of the GenBank output file, defaults\n"
                  + "         to stdout.\n"
                  + "      -f filename\n"
                  + "         Specifies the name of the FASTA output file, defaults\n"
                  + "         to out.fasta.\n"
                  + "      -a [filename]\n"
                  + "         Prints all details of alignment of the GenBank and FASTA\n"
                  + "         sequences, either to the specified file, or to stdout if no\n"
                  + "         filename is given. Off by default.\n");

   }
   public static LinkedHashMap<String, DNASequence> getFASTASequence(
         String filename) {
      File fastaFile = new File(filename);
      try {
         LinkedHashMap<String, DNASequence> seq = FastaReaderHelper
               .readFastaDNASequence(fastaFile);
         if (seq.size() != 1) {
            System.err.printf("FASTA file does not contain a "
                  + "single sequence: %s%n", filename);
            System.exit(1);
         }
         return seq;
      } catch (Exception e) {
         System.err.printf("Incorrectly formatted FASTA file: %s%n",
               filename);
         System.exit(1);
         return null;
      }

   }

   
   public static void writeJoinFeature(genbankFeature feature,
         DNASequence fastaSeq, PrintWriter output,boolean complement){
      int subEnd = feature.description.indexOf("\n");
      String totalLocation = feature.description.substring(0, subEnd);
      //substring from 5 to subEnd-1 so as to get rid of preceding
      //"join(" and trailing ")"
      String featureSequence;
      int numSequence = 1;
      for(int i=0;i<feature.location.length -1;i+=2){
         output.printf(">%s %s Sequence %d %d..%d", feature.type,
               totalLocation,numSequence,feature.location[i],
               feature.location[i+1]);
         if(complement==true)
            featureSequence =fastaSeq.getSubSequence(feature.location[i], 
                  feature.location[i+1])
            .getInverse().getSequenceAsString();
         else
            featureSequence = fastaSeq.getSubSequence(feature.location[i], 
                  feature.location[i+1])
            .getSequenceAsString();
         for (int j = 0; j < featureSequence.length(); ++j) {
            if (j % 60 == 0) output.println();
            output.print(featureSequence.charAt(j));
         }
         output.println();
         ++numSequence;
      }
   }
   
   public static void writeFeatureSequences(List<genbankFeature> features,
         DNASequence fastaSeq, String fastaName){
      PrintWriter fastaOut = null;
      try {
         fastaOut = new PrintWriter(new BufferedWriter(new FileWriter(
               fastaName)));

      } catch (IOException e) {
         System.err.println(e.getMessage());
         System.err.println(e.getStackTrace());
         System.exit(1);
      }

      for (genbankFeature feat : features){
         if(feat.location.length < 2 ||
            !feat.description.substring(0, 
               feat.description.indexOf('\n')).contains(".."))continue;
         if(feat.description.charAt(0) == 'j'){
            writeJoinFeature(feat, fastaSeq, fastaOut,false);
            continue;
         }if(feat.description.charAt(11) == 'j'){
            writeJoinFeature(feat, fastaSeq, fastaOut,true);
            continue;
         }
         String[] qualifiers = feat.description.split("\n");
         if(qualifiers.length > 2) fastaOut.printf(">%s %s %s", 
               feat.type, qualifiers[0],qualifiers[2]);
         else fastaOut.printf(">%s %s", feat.type, qualifiers[0]);
         String featureSequence;
         if (qualifiers[0].charAt(0) == 'c')
            featureSequence = fastaSeq.getSubSequence(feat.location[0], feat.location[1])
            .getInverse().getSequenceAsString();
         else
            featureSequence = fastaSeq.getSubSequence(feat.location[0], feat.location[1])
            .getSequenceAsString();
         for (int i = 0; i < featureSequence.length(); ++i) {
            if (i % 60 == 0)
               fastaOut.println();
            fastaOut.print(featureSequence.charAt(i));
         }
         fastaOut.println();
      }fastaOut.close();
   }


   public static void writeFeatures(List<genbankFeature> features,
         PrintWriter out, int seqlength) {
      out.println("FEATURES             Location/Qualifiers");
      out.printf("     %-16s%s%n", "source", "1.." + seqlength);
      String offset = "                     ";
      int numSpaces = offset.length();
      for (genbankFeature feat : features) {
         String[] qualifiers = feat.description.split("\n");
         out.printf("     %-16s%s%n", feat.type, qualifiers[0]);
         for (int i = 1; i < qualifiers.length; ++i) {
            String[] splitQualifier = qualifiers[i].split("=");
            // Printing type of qualifier first
            int lineLen = numSpaces + splitQualifier[0].length() + 3;
            out.printf("%s/%s=\"", offset, splitQualifier[0]);
            for (int j = 0; j <= splitQualifier[1].length(); ++j) {
               if (lineLen == 79) {
                  out.println();
                  out.print(offset);
                  lineLen = numSpaces;
               }
               if (j == splitQualifier[1].length()) {
                  out.println("\"");
                  break;
               } else
                  out.print(splitQualifier[1].charAt(j));
               ++lineLen;
            }
         }
      }
   }


   public static genbankFeature alignFeature(
         SequencePair<DNASequence, NucleotideCompound> alignment,
         TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> feature) {

      String annotation = feature.getSource();
      //System.out.printf("%n%s => %s%n", feature.getDescription(),annotation);
      int subEnd = annotation.indexOf("\n");
      String location = annotation.substring(0,subEnd);

      String[] nums = location.split("\\D+");
      String[]nonNums = location.split("\\d+");
      boolean numStart = Character.isDigit(location.charAt(0));
      if(numStart) nonNums = Arrays.copyOfRange(nonNums, 1, nonNums.length);
      else nums = Arrays.copyOfRange(nums, 1, nums.length);
      int[] fastaIndices= new int[nums.length];
      int[] gbIndices = new int[nums.length];
      String newLocation="";
      int i = 0;
      for(;i<nums.length;++i){
         try {
            gbIndices[i] = Integer.parseInt(nums[i]);
            fastaIndices[i] = alignment.getIndexInQueryForTargetAt(gbIndices[i]);
         } catch (NumberFormatException e) {
            System.err.printf("Invalid sequence location: %s",
                  annotation.substring(0, subEnd));
            return null;
         }
      }
      //Check for poorly aligned features, if good alignemnt,
      //then count number of gaps introduced in both the genbank
      //and the fasta sequence in this region of the alignment.
      int fastaStart=0;
      int fastaEnd=0;
      int gbStart=0;
      int gbEnd =0;
      int numGaps = 0;
      for(i=0;i<fastaIndices.length -1;i+=2){
         fastaStart=fastaIndices[i];
         fastaEnd=fastaIndices[i+1];
         gbStart=gbIndices[i];
         gbEnd=gbIndices[i+1];
         if (fastaEnd - fastaStart >  1.05 * (gbEnd - gbStart)
                    || fastaEnd-fastaStart < 0.95*(gbEnd-gbStart)) {
                    System.err.println("BAD FEATURE ALIGNMENT");
                    System.err.println("Position in GenBank: "
                            +gbStart+".."+gbEnd
                            +"\nPosition in FASTA: " +fastaStart
                            +".."+fastaEnd);
                    return null;
         }else{
            int j = 0;
            for(NucleotideCompound C:alignment.getQuery()){
               ++j;
               if(j < fastaStart) continue;
               if(j > fastaEnd) break;
               if(C.toString().equals("-")) ++numGaps;

            }j=0;
            for(NucleotideCompound C:alignment.getTarget()){
               ++j;
               if(j < gbStart) continue;
               if(j > gbEnd) break;
               if(C.toString().equals("-")) ++numGaps;
            }
         }
      }
      int max = fastaIndices.length > nonNums.length ? fastaIndices.length : nonNums.length;
      if(numStart){
         for(i =0; i < max;++i){
            if(i<fastaIndices.length)newLocation+=fastaIndices[i];
            if(i< nonNums.length)newLocation += nonNums[i];
         }
      }
      else{
         for(i =0; i < max;++i){
            if(i < nonNums.length)newLocation += nonNums[i];
            if(i<fastaIndices.length)newLocation += fastaIndices[i];
         }
      }
      // New annotation is just the new coordinates concatenated to
      // the old annotation
      String newAnnotation = newLocation;
      if(feature.getDescription().equals("CDS")){
         int index = annotation.indexOf("translation=");
         String CDS;
         if(location.contains("complement"))
            CDS =alignment.getQuery().getOriginalSequence().getSubSequence(fastaStart, fastaEnd).getInverse().getSequenceAsString();
         else CDS =alignment.getQuery().getOriginalSequence().getSubSequence(fastaStart, fastaEnd).getSequenceAsString();
         String newTranslation = new DNASequence(CDS).getRNASequence().getProteinSequence().getSequenceAsString();
         String warning = "";
         if(newTranslation.charAt(0) != 'M') warning+="Inferred CDS does not begin with ATG start codon";
         if(newTranslation.contains("*")) warning+=" Inferred CDS contains premature stop codon.";
         if(!warning.isEmpty())warning = "\nWarning=" + warning.trim();
         newAnnotation += warning+"\nGaps in alignment="
               + numGaps + annotation.substring(subEnd,index)
               +"translation="+newTranslation+'\n'+"evidence=not_experimental\n";
      }else
      newAnnotation = newLocation
               + "\n" + "Gaps in alignment="
               + numGaps + annotation.substring(subEnd);
      genbankFeature genFeat = new genbankFeature(
            feature.getDescription(), newAnnotation, 
            fastaIndices);
      return genFeat;
   }

   public static LinkedList<genbankFeature> alignFeatures(
         GenbankProxySequenceReader<NucleotideCompound> genbankDNAReader,
         SequencePair<DNASequence, NucleotideCompound> alignment,
         LinkedList<String> desiredFeatures) {
      LinkedList<genbankFeature> outputFeatures = new LinkedList<genbankFeature>();
      for (Entry<String, List<TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>>> entry : genbankDNAReader
            .getSequenceParser().mapFeature.entrySet()) {
         if (!desiredFeatures.isEmpty())
            if (!desiredFeatures.contains(entry.getKey()))
               continue;
         for (TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound> feature : entry
               .getValue()) {
            if (feature.getDescription().equals("source"))
               continue;
            genbankFeature feat = alignFeature(alignment,
                  feature);
            if (feat != null) {
               outputFeatures.add(feat);
            }
         }
      }
      return outputFeatures;
   }


   public static SequencePair<DNASequence,NucleotideCompound> getAlignment(
         DNASequence fastaSeq, DNASequence gbSeq, String alignmentFile){
      boolean badAlign=false;
      SequencePair<DNASequence,NucleotideCompound> psa = null;
      PrintWriter alignmentOut = null;
      if(alignmentFile != null){
         if(!alignmentFile.isEmpty()){
            try{
               alignmentOut = new PrintWriter(
                     new BufferedWriter(new FileWriter(alignmentFile)));
            }catch(IOException e){
               System.err.println(e.getMessage());
               System.err.println(e.getStackTrace());
               System.exit(1);
            }
         }else alignmentOut = new PrintWriter(System.out);
         
      }
      int wrapFraction=8;
      int attempts=wrapFraction;
      do{
         if(--attempts == 0) break;
         psa =Alignments.getPairwiseAlignment(
                     fastaSeq,gbSeq,
                     Alignments.PairwiseSequenceAlignerType.GLOBAL,
                     gapP,matrix);
         if(alignmentOut!=null) alignmentOut.println(psa.toString(100));
         String fastaString = fastaSeq.getSequenceAsString();
         boolean rewrap=false;
         boolean gapFasta = true;
         int numNucs=0;
         int numGaps=0;
         boolean foundGap = false;
         int i =1;
         for(;i <= psa.getTarget().getLength();++i){
            if(psa.getTarget().getCompoundAt(i).toString()
                  .equals("-")){
               ++numGaps;
               if(!foundGap) foundGap = true;
            }
            else{
               if(foundGap) break;
               ++numNucs;
            }
         }//If wrapping looked good from beginning of the FASTA sequence
         //then check if it looks good from the end of the FASTA sequence
         //for complete check for wrapping issues.
         if(numNucs > 10 || numGaps < 10){
            gapFasta=false;
            numGaps=numNucs=0;
            foundGap=false;
            for(i=psa.getTarget().getLength();i >0 ;--i){
               if(psa.getTarget().getCompoundAt(i).toString()
                     .equals("-")){
                  ++numGaps;
                  if(!foundGap) foundGap = true;
               }
               else{
                  if(foundGap) break;
                  ++numNucs;
               }
            }
         }if(numNucs <= 10 && numGaps >= 10){
            rewrap = true;
            int index;
            if(!gapFasta)
               index = fastaSeq.getLength()-(numNucs+numGaps);
            else index = numGaps;
            fastaString = fastaString.substring(index)
                  +fastaString.substring(0, index);
            if(alignmentOut!=null)
               alignmentOut.println("Rewrapping FASTA sequence from "
                  + "base " + index);
         }
         if(rewrap){
            DNASequence rewrapSeq = new DNASequence(fastaString, DNACompoundSet.getDNACompoundSet());

            psa =
                  Alignments.getPairwiseAlignment(
                        rewrapSeq,gbSeq,
                        Alignments.PairwiseSequenceAlignerType.GLOBAL,
                        gapP,matrix);
            if(alignmentOut!=null) alignmentOut.println(psa.toString(100));
         }if(psa.getLength() > 1.05*
               (fastaSeq.getLength()>gbSeq.getLength()
                     ? fastaSeq.getLength() : gbSeq.getLength())){
            if(alignmentOut!=null)
               alignmentOut.println("Bad alignment, attempting arbitrary rewrap");
            badAlign=true;
            fastaString  = fastaString.substring(fastaSeq.getLength()/8)
                  +fastaString.substring(0, fastaSeq.getLength()/8);
            fastaSeq=new DNASequence(fastaString, DNACompoundSet.getDNACompoundSet());
         }else badAlign=false;
      }while(badAlign);
      if(alignmentOut != null)alignmentOut.flush();
      //If alignmentOut is not empty, then it's a PrintWriter wrapped around something
      //other than System.out, and hence needs to be closed.
       if(alignmentFile != null && !alignmentFile.isEmpty()) alignmentOut.close();
      return psa;
   }

   public static void main(String[] args){
      String usage = "Usage: mtAnnotate [options] GenBankFile FASTAFile";
      String outGB = null;
      String outFasta = "out.fasta";
      String alignmentFile=null;
      if (args.length < 2) {
         if(args[0].equals("--help")){
            printHelp();
            System.exit(0);
         }
         System.err.println(usage);
         System.err.println("Use option --help for help message.");
         System.exit(1);
      }
      LinkedList<String> wanted_features = new LinkedList<String>(); 
      int argi = 0;
      for(; argi < args.length - 2; ++argi){
         if (args[argi].charAt(0) == '-'){
            char opt = args[argi].charAt(1);
            switch(opt){
            case 't': 
               int j = 1;
               for(;argi+j<args.length-2;++j){
                  if(args[j+argi].charAt(0) =='-')break;
                  wanted_features.add(args[j+argi]);
                  //If j =1 when exiting for loop, then loop only ran through once
               }if(j==1){
                  System.err.println("-f option requires "
                        +"at least one feature type to be "
                        +"specified.");
                  System.err.println(usage);
                  System.exit(1);
               }else argi+=(j-1);
               break;
            case 'g':
               //Check if -g was given argument
               if(argi==args.length-3||args[argi+1].charAt(0)=='-'){
                  System.err.println("-g option requires "
                        +"a filename argument.");
                  System.err.println(usage);
                  System.exit(1);
               }else outGB=args[++argi];
               break;
            case 'f':
               //Check if -f was given argument
               if(argi==args.length-3||args[argi+1].charAt(0)=='-'){
                  System.err.println("-f option requires "
                        +"a filename argument.");
                  System.err.println(usage);
                  System.exit(1);
               }else outFasta=args[++argi];
               break;
            case 'a':
               if(argi==args.length-3||args[argi+1].charAt(0)=='-'){
                  alignmentFile="";
               }else alignmentFile=args[++argi];
               break;
            default:
               System.err.println("Invalid option: " + opt);
            }


         }
      }
      String accession = args[argi];
      //Initialize FASTA file
      String fastaName = args[argi+1];
      LinkedHashMap<String,DNASequence> seq = getFASTASequence(fastaName);
      String fastaAccession ="";
      for (Entry<String, DNASequence> entry: seq.entrySet()){
         fastaAccession = entry.getKey();
      }
      DNASequence fastaSeq= seq.get(fastaAccession);
      File file = new File(accession);
      GenbankProxySequenceReader<NucleotideCompound> genbankDNAReader = null;
      try{
         genbankDNAReader = new GenbankProxySequenceReader
               <NucleotideCompound>(file, accession, DNACompoundSet.getDNACompoundSet());
         
      }catch(Exception e){
         System.err.println(e.getMessage());
         System.err.println(e.getStackTrace());
         System.exit(1);
      }/*
      int numFeats = 0;
      for (Entry<String, List<TextFeature<AbstractSequence<NucleotideCompound>, NucleotideCompound>>> entry : genbankDNAReader
            .getSequenceParser().mapFeature.entrySet()) numFeats+=entry.getValue().size();
      System.out.println("NUMBER FEATURES: "+numFeats);*/
      DNASequence dnaSequence = new DNASequence(genbankDNAReader);
      genbankDNAReader.getHeaderParser().parseHeader(genbankDNAReader.getHeader(), dnaSequence);

      SequencePair<DNASequence,NucleotideCompound> psa = getAlignment(fastaSeq, dnaSequence, alignmentFile);
      fastaSeq = psa.getQuery().getOriginalSequence();
      LinkedList<genbankFeature> outputFeatures = alignFeatures(genbankDNAReader,psa,wanted_features);
      Collections.sort(outputFeatures);
      String[] oldHdr = genbankDNAReader.getHeader().split("\\s+");
      String offset = "            ";
      PrintWriter output = null;
      if(outGB == null) output = new PrintWriter(System.out);
      else {
         try{
            output = new PrintWriter(new BufferedWriter(new FileWriter(outGB)));

         }catch(IOException e){
            System.err.println(e.getMessage());
            System.err.println(e.getStackTrace());
            System.exit(1);
         }
      }
      output.printf("LOCUS       UNDEFINED %18d bp    %s     %s %s %s%n",
            fastaSeq.getLength(),oldHdr[3],oldHdr[4],oldHdr[5],getDate());
      output.printf("DEFINITION  Annotated FASTA sequence from FASTA sequence:%n"
            +offset+"%s%n"
            +offset+"From FASTA file:%n"
            +offset+fastaName +"%n"
            +offset+"Annotations inferred from alignment of features in GenBank sequence:%n"
            +offset+"%s%n", fastaAccession,accession);
      writeFeatures(outputFeatures, output, fastaSeq.getLength());
      writeFeatureSequences(outputFeatures,fastaSeq,outFasta);
      String data = fastaSeq.getSequenceAsString().toLowerCase();
      int seq_len = data.length();
      output.println("ORIGIN");
      int lineLength=60;
      int SEQUENCE_INDENT = 9;
      for (int line_number = 0; line_number < seq_len; line_number += lineLength) {
         output.print(StringManipulationHelper.padLeft(
               Integer.toString(line_number + 1), SEQUENCE_INDENT));
         for (int words = line_number; words < Math.min(line_number
               + lineLength, seq_len); words += 10) {
            if ((words + 10) > data.length()) {
               output.print((" " + data.substring(words)));
            } else {
               output.print((" " + data.substring(words, words + 10)));
            }
         }
         output.println();
      }
      output.println("//");
      output.close();
   }

}
