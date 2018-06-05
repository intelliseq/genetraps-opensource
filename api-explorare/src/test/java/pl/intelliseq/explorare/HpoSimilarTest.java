package pl.intelliseq.explorare;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.junit.Test;

import pl.intelliseq.explorare.model.hpo.HpoOboParser;
import pl.intelliseq.explorare.model.hpo.HpoTerm;
import pl.intelliseq.explorare.model.hpo.HpoTree;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

public class HpoSimilarTest {

	@Test
	public void parserTest() throws IOException {
		HpoTree hpoTree = new HpoOboParser().getHpoTree();
		assertTrue(hpoTree.getTerms().size() > 10000);
		
		System.out.println(hpoTree.getHpoTermById("HP:0002647").isSimilarTo("aorta dissection blabla"));
		
				
	}

	
	
}
