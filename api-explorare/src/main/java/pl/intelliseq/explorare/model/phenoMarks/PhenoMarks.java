package pl.intelliseq.explorare.model.phenoMarks;

import java.util.HashSet;
import java.util.Set;

import com.fasterxml.jackson.annotation.JsonView;

import pl.intelliseq.explorare.model.hpo.HpoTerm;
import pl.intelliseq.explorare.utils.json.Views;

/**
 * @author marpiech
 */
public class PhenoMarks {

	/*** text to parse ***/
	private String rawText;
	
	/*** text with html tags ***/
	@JsonView(Views.Rest.class)
	private StringBuilder markedText = new StringBuilder();
	
	/*** set of identified hpo terms ***/
	@JsonView(Views.Rest.class)
	private Set<HpoTerm> hpoTerms = new HashSet<HpoTerm>();
	
	
	/**
	 * @param rawText text to parse
	 */
	public PhenoMarks(String rawText) {
		this.rawText = rawText;
	}

	public String getRawText() {
		return rawText;
	}

	public String getMarkedText() {
		return markedText.toString();
	}

	public void addMarkedText(String markedText) {
		this.markedText.append(markedText);
	}

	public Set<HpoTerm> getHpoTerms() {
		return hpoTerms;
	}

	public void addHpoTerm(HpoTerm hpoTerm) {
		this.hpoTerms.add(hpoTerm);
	}
	// this is workaround of nltk which does no tag last word.
	// We are adding " end" before processing end remove after 
	public void removeLastWordFromMarkedText() {
		markedText.delete(markedText.lastIndexOf("end") - 1, markedText.length());
		rawText = rawText.substring(0, rawText.lastIndexOf("end") - 1);
	}
	
	
}
