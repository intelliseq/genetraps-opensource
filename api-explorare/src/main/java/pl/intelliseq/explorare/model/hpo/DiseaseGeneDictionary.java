package pl.intelliseq.explorare.model.hpo;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Stream;
import java.util.stream.Collectors;

public class DiseaseGeneDictionary {

	private static String PATH_TO_DISEASES_FILE = "src/main/resources/hpo/phenotype_annotation.tab";
	private static String PATH_TO_GENES_FILE = "src/main/resources/hpo/diseases_to_genes.txt";
	
	private Map<String, Disease> diseases = new TreeMap<String, Disease>();
	
	public DiseaseGeneDictionary() {
		this.parseDiseasesFile();
		this.parseGenesFile();
	}
	
	public boolean add(Disease disease) {
		if (diseases.containsKey(disease.getId())) return false;
		this.diseases.put(disease.getId(), disease);
		return true; }
	
	public Disease getDiseaseById(String id) {
		return this.diseases.get(id); }
	
	public Disease getDiseaseByName(String diseaseName) {
		return this.diseases.values().stream()
				.filter(disease -> disease.getPrimaryName().equals(diseaseName))
				.findFirst()
	            .get();
	}
	
	public Set<String> getDiseases() {
		return this.diseases.entrySet().stream()
				.map(Map.Entry::getValue)
				.map(Disease::getPrimaryName)
				.collect(Collectors.toSet());
	}
	
	private void parseDiseasesFile() {
		try (Stream<String> stream = Files.lines(Paths.get(PATH_TO_DISEASES_FILE))) {

			stream. forEach(s -> 
				{
					String[] elements = s.split("\t");
					String diseaseDatabase = elements[0];
					String diseaseDatabaseId = elements[1];
					String diseaseId = diseaseDatabase + ":" + diseaseDatabaseId;
					Disease disease = new Disease(diseaseId);
					if(this.add(disease)) {
						for(String diseaseName : elements[2].split(";")) {
							diseaseName = diseaseName.replaceAll("^?[0-9]{6} ","").replaceAll("[\\*|#|%|\\+]", "");
							if(!diseaseName.contains("MOVED TO")) {
								disease.setPrimaryName(diseaseName);
								break;
							}
						}
						disease.setSynonyms(
								new HashSet<String>(Arrays.asList(elements[2].split(";"))));
						if ((disease.getPrimaryName() == null))
							this.diseases.remove(diseaseId);
							//System.out.println(elements[2]);
						/*disease.setPrimaryName(
								elements[2]
									.split(";")[0]
									.replaceAll(".[0-9]{6} ","")
									);*/
									//.replaceAll("#[0-9]+ ", "")
									//.replaceAll("\\+[0-9]+ ", "")
									//.replaceAll("%[0-9]+ ", ""));
						
					}
				});
			
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	private void parseGenesFile() {
		try (Stream<String> stream = Files.lines(Paths.get(PATH_TO_GENES_FILE))) {

			stream. forEach(s -> 
				{
					String[] elements = s.split("\t");
					if (elements.length > 1) {
						//long startTime = System.currentTimeMillis();
						Disease disease = this.getDiseaseById(elements[0]);
						//long secondTime = System.currentTimeMillis();
						if (disease != null) disease.addGene(elements[2]);
						//long endTime = System.currentTimeMillis();
						//System.out.println("Times " + (secondTime - startTime) + " : " + (endTime - secondTime));
					}
				});
			
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}


	
}
