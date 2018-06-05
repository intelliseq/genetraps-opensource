package pl.intelliseq.explorare.mvc;

import java.util.List;

import pl.intelliseq.explorare.model.hpo.HpoTerm;

public class PhenotypeQuery {

	public List<HpoTerm> hpoTerms;

	public List<HpoTerm> getHpoTerms() {
		return hpoTerms;
	}

	public void setHpoTerms(List<HpoTerm> hpoTerms) {
		this.hpoTerms = hpoTerms;
	}
	
	
	
}
