package pl.intelliseq.explorare.model.hpo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

public class HpoTree {
	
	HpoTerm root;
	Map <String, HpoTerm> hpoMap = new TreeMap <String, HpoTerm> ();
	private Set <String> cachedDiseases;
	
	public HpoTree() {
		root = new HpoTerm("HP:0000001");
		hpoMap.put(root.getId(), root);
	}
	
	public void add(List<String> lineBuffer) {
		
		HpoTerm term = null;
		for (String line : lineBuffer) {
			
			String[] elements = line.split(": ");
			
			if (elements[0].equals("id"))
				term = this.getOrCreateHpoTermById(elements[1]);

			if (elements[0].equals("name"))
				term.setName(elements[1]);
			
			if (elements[0].equals("synonym"))
				term.addSynonym(elements[1].split("\"")[1]);
			
			if (elements[0].equals("is_a"))
				term.setParent(
						this.getOrCreateHpoTermById(elements[1].split(" ")[0])
						);
				
		} // end parsing line from lineBuffer
	}

	public HpoTerm getOrCreateHpoTermById(String id) {
		return this.getOrCreateHpoTermById(id, true);
	}
	
	public HpoTerm getHpoTermById(String id) {
		return this.getOrCreateHpoTermById(id, false);
	}	
	
	public void calculateWeights() {
		this.getHpoTermById("HP:0000001").setWeight(1d);
	}
	
	private HpoTerm getOrCreateHpoTermById(String id, boolean create) {
		if (hpoMap.containsKey(id)) return hpoMap.get(id);
		
		if (create) {
			HpoTerm newTerm = new HpoTerm(id);
			hpoMap.put(id, newTerm);
			return newTerm;
		}
		
		throw new RuntimeException("No such HPO term " + id);
	}
	
	
	
	public List<HpoTerm> getTerms() {
		return new ArrayList<HpoTerm>(hpoMap.values());
	}

	public List<HpoTerm> getTermsForOboParsing() {
		List <HpoTerm> termsForParsing = new ArrayList<HpoTerm>();
		for (HpoTerm term : hpoMap.values())
			if (term.hasParent("HP:0000118"))
				termsForParsing.add(term);
		return termsForParsing;
	}
	
	public void addGene(String id, String gene) {
		HpoTerm term = this.getHpoTermById(id);
		term.addGene(gene);
	}

	public void addDisease(String id, String disease) {
		HpoTerm term = this.getHpoTermById(id);
		term.addDisease(disease);
	}
	
	public Set<String> getDiseases() {
		if (cachedDiseases == null) {
			//System.out.println("NULL");
			Set <String> diseases = new HashSet <String>();
			for(HpoTerm term : this.getTerms())
				diseases.addAll(term.getDiseases());
			this.cachedDiseases = diseases;
		}
		//System.out.println(this.cachedDiseases.size());
		return cachedDiseases;
	}
	
	public Map <String, Double> getGenes(Set<HpoTerm> hpoTerms) {
		
		Map <String, Double> nonNormalizedGeneScores = new LinkedHashMap <String, Double>();
		
		for (HpoTerm symptom : hpoTerms) {
			//HpoTerm symptom = this.getHpoTermById(symptomId);
			symptom.increaseGeneScoresByWeight(nonNormalizedGeneScores, 1d);
			try {symptom.getParent().increaseGeneScoresByWeight(nonNormalizedGeneScores, 0.5d);} catch (Exception e) {}
			try {symptom.getParent().getParent().increaseGeneScoresByWeight(nonNormalizedGeneScores, 0.25d);} catch (Exception e) {}
		}
		
		Map <String, Double> sortedMap = HpoTree.sortByValue(nonNormalizedGeneScores);
		
		/* normalize sorted map */
		Double maximum = sortedMap.remove("maximum");
		sortedMap.entrySet().forEach(entry -> entry.setValue(entry.getValue() / maximum));
		
		return sortedMap;
		
	}
	
	public Map <String, Double> getDiseases(Set<HpoTerm> hpoTerms) {
		
		Map <String, Double> nonNormalizedDiseaseScores = new LinkedHashMap <String, Double>();
		
		for (HpoTerm symptom : hpoTerms) {
			//HpoTerm symptom = this.getHpoTermById(symptomId);
			symptom.increaseDiseaseScoresByWeight(nonNormalizedDiseaseScores, 1d);
			try {symptom.getParent().increaseDiseaseScoresByWeight(nonNormalizedDiseaseScores, 0.5d);} catch (Exception e) {}
			try {symptom.getParent().getParent().increaseDiseaseScoresByWeight(nonNormalizedDiseaseScores, 0.25d);} catch (Exception e) {}
		}
		
		Map <String, Double> sortedMap = HpoTree.sortByValue(nonNormalizedDiseaseScores);
		
		/* normalize sorted map */
		Double maximum = sortedMap.remove("maximum");
		sortedMap.entrySet().forEach(entry -> entry.setValue(entry.getValue() / maximum));
		
		return sortedMap;
		
	}
	
	private static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
	    return map.entrySet()
	              .stream()
	              .sorted(Map.Entry.comparingByValue(Collections.reverseOrder()))
	              .collect(Collectors.toMap(
	                Map.Entry::getKey, 
	                Map.Entry::getValue, 
	                (e1, e2) -> e1, 
	                LinkedHashMap::new
	              ));
	}

	public Set<HpoTerm> getPhenotypesByDisease(String diseaseName) {
		Set<HpoTerm> hpoTerms = new HashSet<HpoTerm>();
		for(HpoTerm term : this.getTerms()) {
			if (term.causeDisease(diseaseName))
				hpoTerms.add(term);
		}
		return hpoTerms;
	}

	


	
}
