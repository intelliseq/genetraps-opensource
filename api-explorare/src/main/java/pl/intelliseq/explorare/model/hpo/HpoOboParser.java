package pl.intelliseq.explorare.model.hpo;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class HpoOboParser {

	private static String PATH_TO_OBO_FILE = "src/main/resources/hpo/hp.obo";
	private static String PATH_TO_GENES_FILE = "src/main/resources/hpo/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt";
	private static String PATH_TO_DISEASES_FILE = "src/main/resources/hpo/phenotype_annotation.tab";
	private HpoTree hpoTree = new HpoTree();
	
	public HpoOboParser () {
		
		try (Stream<String> stream = Files.lines(Paths.get(PATH_TO_OBO_FILE))) {

			HpoCollector collector = new HpoCollector(hpoTree);
			stream. forEach(s -> collector.parse(s));
			
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		try (Stream<String> stream = Files.lines(Paths.get(PATH_TO_GENES_FILE))) {

			GeneCollector geneCollector = new GeneCollector(hpoTree);
			stream. forEach(s -> geneCollector.parse(s));
			
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		try (Stream<String> stream = Files.lines(Paths.get(PATH_TO_DISEASES_FILE))) {

			DiseaseCollector diseaseCollector = new DiseaseCollector(hpoTree);
			stream. forEach(s -> diseaseCollector.parse(s));
			
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		hpoTree.calculateWeights();
		
	}

	public List<HpoTerm> getTerms() {
		return hpoTree.getTerms();
	}
	
	public List<HpoTerm> getTermsForOboParsing() {
		return hpoTree.getTermsForOboParsing();
		
	}
	
	public HpoTree getHpoTree() {
		return hpoTree;
	}
	
}

class HpoCollector {
	
	private boolean isReading = false;
	private HpoTree hpoTree;
	private List <String> lineBuffer = new ArrayList <String> ();
	
	public HpoCollector(HpoTree hpoTree) {
		this.hpoTree = hpoTree;
	}

	public void parse(String line) {
		if (isReading) addToTerm(line);
		else if (line.equals("[Term]")) {
			isReading = true;
			lineBuffer = new ArrayList <String> ();
		}
	}
	
	void addToTerm(String line) {
		if (line.isEmpty()) {
			isReading = false;
			hpoTree.add(lineBuffer);
			return;
		} else {
			lineBuffer.add(line);
		}
	}

}

class GeneCollector {
	private HpoTree hpoTree;
	
	public GeneCollector(HpoTree hpoTree) {
		this.hpoTree = hpoTree;
	}

	public void parse(String line) {
		if (line.startsWith("#")) return;
		String[] elements = line.split("\t");
		hpoTree.addGene(elements[0], elements[3]);
	}
}

class DiseaseCollector {
	private HpoTree hpoTree;
	
	public DiseaseCollector(HpoTree hpoTree) {
		this.hpoTree = hpoTree;
	}

	public void parse(String line) {
		//if (line.startsWith("#")) return;
		String[] elements = line.split("\t");
		String diseaseDatabase = elements[0];
		String diseaseId = elements[1];
		String diseaseName = diseaseDatabase + ":" + diseaseId;
		//String diseaseName = elements[2].split(";")[0].replaceAll("#[0-9]+ ", "");
		//if (!diseaseName.equals("MOVED TO"))
			hpoTree.addDisease(elements[4], diseaseName);
	}
}

