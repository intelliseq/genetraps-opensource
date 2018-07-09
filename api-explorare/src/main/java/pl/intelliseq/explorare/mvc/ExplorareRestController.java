package pl.intelliseq.explorare.mvc;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.RequestParam;
import org.springframework.web.bind.annotation.RestController;

import com.fasterxml.jackson.annotation.JsonView;

import pl.intelliseq.explorare.model.Results;
import pl.intelliseq.explorare.model.hpo.DiseaseGeneDictionary;
import pl.intelliseq.explorare.model.hpo.DiseaseResult;
import pl.intelliseq.explorare.model.hpo.GeneResult;
import pl.intelliseq.explorare.model.hpo.HpoTerm;
import pl.intelliseq.explorare.model.hpo.HpoTree;
import pl.intelliseq.explorare.model.hpo.PhenotypeResult;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarks;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;
import pl.intelliseq.explorare.mvc.response.HelloRestResponse;
import pl.intelliseq.explorare.utils.json.Views;

@RestController
public class ExplorareRestController {

    @Autowired
    PhenoMarksParser phenoMarksParser;
    
    @Autowired
    HpoTree hpoTree;
    
    @Autowired
    DiseaseGeneDictionary diseaseGeneDictionary;
    
    @CrossOrigin()
    @RequestMapping(path = "/hello", method = RequestMethod.GET, produces = "application/json")
    public HelloRestResponse greeting() { return new HelloRestResponse("explorare service is alive"); }
    
    @CrossOrigin()
    @JsonView(Views.Rest.class)
    @RequestMapping(path = "/parse-text", method = RequestMethod.POST)
    public PhenoMarks parseText(@RequestBody ParseTextQuery query) {
    	
        PhenoMarks phenoMarks = phenoMarksParser.tagInput(query.getQuery());
        return phenoMarks;
                            
    }
  
    static class GetGenesRequest { 
    	Set<HpoTerm> hpoTerms;
    	Double threshold = 0.3d;
		public void setHpoTerms(Set<HpoTerm> hpoTerms) {this.hpoTerms = hpoTerms;}
		public void setThreshold(Double threshold) {this.threshold = threshold;}
    	
    }
    
    @CrossOrigin()
    @RequestMapping(path = "/get-diseases", method = RequestMethod.POST, consumes = {"application/json"})
    public Set<DiseaseResult> getDiseases(@RequestBody GetGenesRequest request) {
    	Set<HpoTerm> hpoTerms = this.getParsedHpoTerms(request.hpoTerms);
    	
    	//System.out.println("Diseases " + this.hpoTree.getDiseases(hpoTerms).size());
    	
    	Map<String, Double> result = this.hpoTree.getDiseases(hpoTerms)
    		.entrySet()
    		.stream()
    		.filter(map -> map.getValue() > request.threshold)
    		.collect(Collectors.toMap(p -> p.getKey(), p -> p.getValue()));
    	
    	//System.out.println("Result " + result.size());
    	
    	// sorting
    	result = result.entrySet().stream()
        .sorted(Map.Entry.<String, Double>comparingByValue().reversed())
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                (e1, e2) -> e1, LinkedHashMap::new));
    	
    	//System.out.println("Result " + result.size());
    	
    	Set<DiseaseResult> finalResult = new TreeSet<DiseaseResult>();
    	for(Entry<String, Double> entry : result.entrySet()) {
    		//System.out.println(entry.getKey());
    		if (this.diseaseGeneDictionary.getDiseaseById(entry.getKey()) != null)
    		finalResult.add(
    				new DiseaseResult(
    						this.diseaseGeneDictionary
    						.getDiseaseById(entry.getKey())
    						.getPrimaryName(),
    						100 * entry.getValue()
    						)
    				);
    	}
    	//System.out.println("Final Result " + finalResult.size());
    	
    	return finalResult;
                            
    }
    
    @CrossOrigin()
    @RequestMapping(path = "/get-genes", method = RequestMethod.POST, consumes = {"application/json"})
    public Set<GeneResult> getGenes(@RequestBody GetGenesRequest request) {
    	Set<HpoTerm> hpoTerms = this.getParsedHpoTerms(request.hpoTerms);

    	// getting diseases
    	Map<String, Double> result = this.hpoTree.getDiseases(hpoTerms)
    		.entrySet()
    		.stream()
    		.filter(map -> map.getValue() > request.threshold)
    		.collect(Collectors.toMap(p -> p.getKey(), p -> p.getValue()));
    	
    	// sorting
    	result = result.entrySet().stream()
        .sorted(Map.Entry.<String, Double>comparingByValue().reversed())
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                (e1, e2) -> e1, LinkedHashMap::new));
    	
	
    	Set<GeneResult> finalResult = new TreeSet<GeneResult>();
    	for(Entry<String, Double> entry : result.entrySet()) {
    		if (this.diseaseGeneDictionary.getDiseaseById(entry.getKey()) != null)
    			for (String gene : this.diseaseGeneDictionary.getDiseaseById(entry.getKey()).getGenes())
    				finalResult.add(new GeneResult(gene, 100 * entry.getValue()));
    				
    	}
    	//System.out.println("Final Result " + finalResult.size());
    	
    	return finalResult;
                            
    }
    
    private Set <HpoTerm> getParsedHpoTerms(Set <HpoTerm> simpleHpoTerms) {
    	Set <HpoTerm> parsedHpoTerms = new HashSet<HpoTerm>();
    	for (HpoTerm simpleHpoTerm : simpleHpoTerms)
    		parsedHpoTerms.add(this.hpoTree.getHpoTermById(simpleHpoTerm.getId()));
    	return parsedHpoTerms;
    }
    
    @CrossOrigin()
    @JsonView(Views.Rest.class)
    @RequestMapping(path = "/get-phenotypes-by-disease", method = RequestMethod.GET)
    public Results getPhenotypesByDisease(@RequestParam String disease) {
    	
    	String diseaseId = this.diseaseGeneDictionary.getDiseaseByName(disease).getId();
    	
    	Results results = new Results();
    	
    	results.addResult(
    			this.makeSetSpecific(
    					this.hpoTree.getPhenotypesByDisease(diseaseId)
    			));
    	
        return results;
                            
    }
    
    private Set<HpoTerm> makeSetSpecific(Set<HpoTerm> hpoTerms) {
    	Set<HpoTerm> out = new HashSet<HpoTerm>();
    	for(HpoTerm term : hpoTerms) {
    		boolean hasChildren = false;
    		for(HpoTerm childTerm : hpoTerms) {
    			if(childTerm.isChildOf(term))
    				hasChildren = true;
    		}
    		if(!hasChildren)
    			out.add(term);
    	}
    	return out;
    }
    
    @CrossOrigin()
    @JsonView(Views.Rest.class)
    @RequestMapping(path = "/disease-autocomplete", method = RequestMethod.GET)
    public Results diseasesAutocomplete(
    		@RequestParam String firstLetters,
    		@RequestParam Integer resultsCount) {
    	
    	Results results = new Results();
    	
    	for(String disease: this.diseaseGeneDictionary.getDiseases()) {
    		if(disease.toLowerCase().startsWith(firstLetters.toLowerCase())) {
    			results.addResult(disease);
    			if(results.sizeGreaterOrEqualTo(resultsCount))
    				return results;
    		}
    	}
    	
    	for(String disease: this.diseaseGeneDictionary.getDiseases()) {
    		if(disease.toLowerCase().contains(firstLetters.toLowerCase())) {
    			results.addResult(disease);
    			if(results.sizeGreaterOrEqualTo(resultsCount))
    				return results;
    		}
    	}
    	
        return results;//phenoMarks;
                            
    }
    
    @CrossOrigin()
    @JsonView(Views.Rest.class)
    @RequestMapping(path = "/get-phenotype-by-id", method = RequestMethod.GET)
    public PhenotypeResult getPhenotypeById(@RequestParam String phenotypeid) {
    	return new PhenotypeResult(this.hpoTree.getHpoTermById(phenotypeid)); }
    
    @CrossOrigin()
    @JsonView(Views.Rest.class)
    @RequestMapping(path = "/phenotype-autocomplete", method = RequestMethod.GET)
    public Results phenotypeAutocomplete(
    		@RequestParam String firstLetters,
    		@RequestParam Integer resultsCount) {
    	
    	Results results = new Results();
    	
    	for(HpoTerm hpoTerm : this.hpoTree.getTerms()) {
    		if(hpoTerm.getName().toLowerCase().startsWith(firstLetters.toLowerCase())) {
    			results.addResult(hpoTerm);
    			if(results.sizeGreaterOrEqualTo(resultsCount))
    				return results;
    		}
    	}
    	
    	for(HpoTerm hpoTerm : this.hpoTree.getTerms()) {
    		if(hpoTerm.getName().toLowerCase().contains(firstLetters.toLowerCase())) {
    			results.addResult(hpoTerm);
    			if(results.sizeGreaterOrEqualTo(resultsCount))
    				return results;
    		}
    	}
    	
        return results;
                            
    }
    
}