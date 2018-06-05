package pl.intelliseq.explorare.model.hpo;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.DeflaterInputStream;

import com.fasterxml.jackson.annotation.JsonView;

import info.debatty.java.stringsimilarity.Damerau;
import pl.intelliseq.explorare.utils.json.Views;

/**
 * @author marpiech
 *
 */



public class HpoTerm {

	@JsonView(Views.Rest.class)
	String id;
	@JsonView(Views.Rest.class)
	String name;
	List<String> synonyms;
    List<HpoTerm> parents;
	List <HpoTerm> children;
	Set <String> genes;
	Set <String> diseases;
	List <ObsoleteHpoTerm> replacementOf;


	/*** what is weight? Baby don't hurt me. No more.
	 		dear future developer sorry I don't remember ***/

	Double weight;
	
	public HpoTerm() {
		this.id = "";
		this.name = "";
		this.synonyms = new ArrayList <String> ();
		this.children = new ArrayList <HpoTerm> ();
        this.parents = new ArrayList <HpoTerm> ();
		this.genes = new HashSet<String> ();
		this.diseases = new HashSet<String> ();
		this.replacementOf = new ArrayList<>();

	}
	
	public HpoTerm(String id) {
		
		if (!id.matches("HP:[0-9]{7}"))
			throw new RuntimeException("Wrong id format: " + id);
		
		this.id = id;
		this.synonyms = new ArrayList <String> ();
		this.children = new ArrayList <HpoTerm> ();
        this.parents = new ArrayList <HpoTerm> ();
		this.genes = new HashSet<String> ();
		this.diseases = new HashSet<String> ();
        this.replacementOf = new ArrayList<>();

	}

	public void setId(String id) {
		this.id = id;
	}
	
	public void setName(String name) {
		this.name = name;
	}

	public void addParent(HpoTerm parent) {
        parent.addChild(this);
        this.parents.add(parent);
	}
	
	public void addSynonym(String synonym) {
		this.synonyms.add(synonym);
	}

	public void addChild(HpoTerm child) {
		this.children.add(child);
	}
	
	public void addGene(String gene) {
		this.genes.add(gene);
        if (! this.parents.isEmpty()) {
            for (HpoTerm parent : this.parents) {
                parent.addGene(gene);
            }
        }
	}

	public void addDisease(String disease) {
		this.diseases.add(disease);
        if (! this.parents.isEmpty()) {
            for (HpoTerm parent : this.parents){
                parent.addDisease(disease);
            }
        }
	}
	
	public void setWeight(double weight) {
		this.weight = weight;
		for(HpoTerm child : children) {
			child.setWeight(weight + 1d);
		}
	}
	
	public String getId() {
		return id;
	}

	public String getName() {
		return name;
	}

    public Set <String> getGenes() {
		return genes;
	}
	
	public Set <String> getDiseases() {
		return diseases;
	}
	
	public boolean isEqualTo(String phrase) {
		return isSimilarTo(phrase, 0);
	}
	
	public boolean isSimilarTo(String phrase) {
		return isSimilarTo(phrase, 5);
	}

	public boolean isSimilarTo(String phrase, Integer allowedBasesPerMismatch) {
	
		Integer allowedMismatches = 0;
		if (allowedBasesPerMismatch > 0)
		allowedMismatches = (int) (((float) phrase.length())/ ((float) allowedBasesPerMismatch));
		//System.out.println(phrase + allowedMismatches);
		
		phrase = phrase.replaceAll("[^a-zA-Z ]", "").toLowerCase();
		try {
			if( phrase.equals(name.toLowerCase())) return true;
		} catch (Exception e) {
			System.out.println(phrase);
			System.out.println(this.id);
			System.out.println(this.name);
		}
		
		if(allowedMismatches > 0 && phrase.startsWith(name.toLowerCase().substring(0, 3))) {
			Damerau distanceCalculator = new Damerau();
			double distance = distanceCalculator.distance(phrase, name.toLowerCase());
			if (distance <= allowedMismatches) return true;
		}

		if (synonyms != null && synonyms.size() > 0)
			for(String synonym : synonyms)
				if( phrase.equals(synonym.toLowerCase())) return true;
		
		return false;
			
	}

	public List<HpoTerm> getChildren() {
		return children;
	}

	public double getGeneWeight() {
		if (this.genes.isEmpty()) return 0d;
		return 1d / (double) this.genes.size();
	}

	public List <HpoTerm> getParents() {
		return this.parents;
	}

    public List<ObsoleteHpoTerm> getReplacementOf() {
        return replacementOf;
    }

    public void addReplacementOf(String id, String name, List <String> synonyms) {
        ObsoleteHpoTerm obsoleteHpoTerm = new ObsoleteHpoTerm(id, name, synonyms);
        this.replacementOf.add(obsoleteHpoTerm);
    }


    public void increaseGeneScoresByWeight(Map<String, Double> geneScores, double weight) {


        if (this.weight != null) {

            double increaseBy = this.weight * weight;

            for (String symbol : this.genes)
                if (geneScores.containsKey(symbol))
                    geneScores.put(symbol, geneScores.get(symbol) + increaseBy);
                else
                    geneScores.put(symbol, increaseBy);

            if (geneScores.containsKey("maximum"))
                geneScores.put("maximum", geneScores.get("maximum") + increaseBy);
            else
                geneScores.put("maximum", increaseBy);
        }

	}

	public void increaseDiseaseScoresByWeight(Map<String, Double> diseaseScores, double weight) {
		
		double increaseBy = this.weight * weight;
		
		for (String disease : this.diseases)
			if (diseaseScores.containsKey(disease))
				diseaseScores.put(disease, diseaseScores.get(disease) + increaseBy);
			else
				diseaseScores.put(disease, increaseBy);
		
		if (diseaseScores.containsKey("maximum"))
			diseaseScores.put("maximum", diseaseScores.get("maximum") + increaseBy);
		else
			diseaseScores.put("maximum", increaseBy);
		
	}

	public boolean causeDisease(String diseaseName) {
		for (String disease : this.diseases)
			if(disease.equals(diseaseName))
				return true;
		return false;
	}

	public boolean isChildOf(HpoTerm term) {
		return checkIfTermIsAmong(term, this.parents);
	}

	private boolean checkIfTermIsAmong(HpoTerm termToCheck, List <HpoTerm> terms) {
	    if (terms.isEmpty()) return false;
	    for (HpoTerm term : terms) {
	        if (term.equals(termToCheck))
	            return true;
	    }
        List <HpoTerm> parentsOfTerms = new ArrayList<>();
        for (HpoTerm term : terms) {
            parentsOfTerms.addAll(term.getParents());
        }
        return checkIfTermIsAmong(termToCheck, parentsOfTerms);
    }


	public boolean hasParent(String parentId) {
		return checkByIdIfTermIsAmong(parentId, this.parents);
	}

    private boolean checkByIdIfTermIsAmong(String idToCheck, List <HpoTerm> terms) {
        if (terms.isEmpty()) return false;
        for (HpoTerm term : terms) {
            if (term.getId().equals(idToCheck))
                return true;
        }
        List <HpoTerm> parentsOfTerms = new ArrayList<>();
        for (HpoTerm term : terms) {
            parentsOfTerms.addAll(term.getParents());
        }
        return checkByIdIfTermIsAmong(idToCheck, parentsOfTerms);
    }

}
