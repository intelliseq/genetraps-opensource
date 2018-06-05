package pl.intelliseq.explorare.model.phenoMarks;

import opennlp.tools.cmdline.postag.POSModelLoader;
import opennlp.tools.postag.POSModel;
import opennlp.tools.postag.POSSample;
import opennlp.tools.postag.POSTaggerME;
import opennlp.tools.tokenize.WhitespaceTokenizer;
import opennlp.tools.util.InputStreamFactory;
import opennlp.tools.util.ObjectStream;
import opennlp.tools.util.PlainTextByLineStream;
import pl.intelliseq.explorare.model.hpo.HpoOboParser;
import pl.intelliseq.explorare.model.hpo.HpoTerm;
import pl.intelliseq.explorare.model.hpo.HpoTree;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import org.apache.commons.io.IOUtils;

import com.sun.scenario.effect.impl.sw.sse.SSEBlend_SRC_OUTPeer;

public class PhenoMarksParser {

	POSTaggerME tagger;
	List <HpoTerm> terms;
	
	public PhenoMarksParser() {
		POSModel model = new POSModelLoader().load(new File("src/main/resources/dict/en-pos-maxent.bin"));
	    this.tagger = new POSTaggerME(model);
	    this.terms = new HpoOboParser().getTermsForOboParsing();
	}
	
	public PhenoMarksParser(HpoTree hpoTree) {
		POSModel model = new POSModelLoader().load(new File("src/main/resources/dict/en-pos-maxent.bin"));
	    this.tagger = new POSTaggerME(model);
	    //this.terms = hpoTree.getTerms();
	    this.terms = hpoTree.getTermsForOboParsing();
	}

	public PhenoMarks tagInput(String input) {
		
		try {
			return this.tagInputWithCheckedException(input);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	private PhenoMarks tagInputWithCheckedException(String input) throws IOException {
		
		
		/*** define charset ***/
		Charset charset = Charset.forName("UTF-8");
		
		/*** create inputStream ***/
	    ObjectStream<String> lineStream =
	            new PlainTextByLineStream(new InputStreamFactory() {
	                @Override
	                public InputStream createInputStream() throws IOException {
	                    return IOUtils.toInputStream(input + " end", charset);
	                }
	            }, charset);
        
        /*** prepare output ***/
		PhenoMarks phenoMarks = new PhenoMarks(input + " end");
		
	    /* reading lines from input */
	    String line;
	    while ((line = lineStream.read()) != null) {

	        String whitespaceTokenizerLine[] = WhitespaceTokenizer.INSTANCE.tokenize(line);
	        String[] tags = tagger.tag(whitespaceTokenizerLine);
	        
	        POSSample sample = new POSSample(whitespaceTokenizerLine, tags);
	        
	        List <String> tagList = Arrays.asList(sample.getTags());
	        List <String> wordList = Arrays.asList(sample.getSentence());
	        
	        int startTagging = 0;
	        for (int iter = 0; iter < tagList.size(); iter++) {
	        	
	        	//System.out.print(wordList.get(iter) + ":");
	        	//System.out.print(tagList.get(iter) + ": ");
	        	//System.out.println();
	        	String tag = tagList.get(iter);
	        	/* copy text */
	        	phenoMarks.addMarkedText(wordList.get(iter) + " ");
	        	
	        	HpoTerm bestMatch = null;
	        	if (tag.matches("NN|NNP|NNS|JJ|VBD|IN|DT")) {
	        		
	        		String phrase = "";
	        		LinkedList <String> words = new LinkedList <String> ();
	        		
	        		/*** phrase maximum lenght is 4 ***/
	        		int PHRASE_MAXIMUM_LENGTH = 4;
	        		if(startTagging < iter - PHRASE_MAXIMUM_LENGTH + 1)
	        			startTagging = iter - PHRASE_MAXIMUM_LENGTH + 1;
	        		
	        		/*** get list of words ***/
	        		for (int back = 0; back <= (iter - startTagging); back++)
	        			words.addFirst(wordList.get(iter - back)
	        					.replaceAll("[^A-Za-z0-9 ]", ""));
	        					//.replaceAll("\\.|\\,|\\(|\\)", ""));

	        		//for(String word : words) System.out.print(word + ", ");
	        		//System.out.println();
	        		
	        		/*** search for phrases ***/
	        		for (int back = 0; back <= (iter - startTagging); back++) {
	        			phrase = wordList.get(iter - back)
	        					.replaceAll("\\.|\\,|\\(|\\)", "")
	        					.concat(" ")
	        					.concat(phrase)
	        					.trim();
	        			/*** exact matches ***/
	        			for (HpoTerm term : terms)
	        				if (term.isEqualTo(phrase)) {
	        					bestMatch = term;
	        				}
	        			/*** partial matches ***/
	        			if (bestMatch == null)
	        			for (HpoTerm term : terms)
	        				if (term.isSimilarTo(phrase)) {
	        					bestMatch = term;
	        				}
	        			/*** check of the ***/
	        			//System.out.println(phrase);
	        		}
	        		
	        		/*** find pairs of words that match ***/
	        		if (bestMatch == null)
	        		if (words.size() > 2) {
	        			//System.out.println(words.getLast() + " " + words.getFirst());
	        			//if (words.size() > 3) words.pollFirst();
	        			/*** exact matches ***/
	        			for (HpoTerm term : terms)
	        				if (term.isEqualTo(words.getFirst() + " " + words.getLast()))
	        					bestMatch = term;
	        			/*** partial matches ***/
	        			if (bestMatch == null)
	        			for (HpoTerm term : terms)
	        				if (term.isSimilarTo(words.getFirst() + " " + words.getLast()))
	        					bestMatch = term;
	        			if (bestMatch == null)
		        			for (HpoTerm term : terms)
		        				if (term.isSimilarTo(words.getLast() + " " + words.getFirst()))
		        					bestMatch = term;
	        			words.pollFirst();
	        			if (bestMatch == null)
		        			for (HpoTerm term : terms)
		        				if (term.isSimilarTo(words.getFirst() + " " + words.getLast()))
		        					bestMatch = term;
		        			
	        		}
	        		
	        	}
	        	if (bestMatch != null) {
	        		phenoMarks.addMarkedText(("<span class=\"badge badge-info\">" + bestMatch.getId() + "::" + bestMatch.getName() + "</span> "));
	        		phenoMarks.addHpoTerm(bestMatch);
	        	}

	        	if (!(tag.matches("NN|NNP|NNS|JJ|VBD|IN|DT") )) startTagging = iter + 1;
	        }

	        
	    }
	    lineStream.close();
	    phenoMarks.removeLastWordFromMarkedText();
	    return phenoMarks;
	    
	}
	
}
