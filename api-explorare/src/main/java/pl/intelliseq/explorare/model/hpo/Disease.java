package pl.intelliseq.explorare.model.hpo;

import java.util.HashSet;
import java.util.Set;

public class Disease {
	private String id;
	private String primaryName;
	private Set<String> synonyms;
	private Set<String> genes = new HashSet<String>();
	
	public Disease(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getPrimaryName() {
		return primaryName;
	}

	public void setPrimaryName(String primaryName) {
		this.primaryName = primaryName;
	}

	public Set<String> getSynonyms() {
		return synonyms;
	}

	public void setSynonyms(Set<String> synonyms) {
		this.synonyms = synonyms;
	}

	public Set<String> getGenes() {
		return genes;
	}

	public void setGenes(Set<String> genes) {
		this.genes = genes;
	}
	
	public boolean addGene(String gene) {
		return this.genes.add(gene);
	}
	
	@Override
	public boolean equals(Object obj) {
	    if (obj == null) {
	        return false;
	    }
	    if (!Disease.class.isAssignableFrom(obj.getClass())) {
	        return false;
	    }
	    final Disease other = (Disease) obj;
	    if (this.id == null) return false;
	    if (other.id == null) return false;
	    if (this.id.equals(other.id)) return true;
	    return false;
	}
	
 }
