package pl.intelliseq.explorare;

import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Test;

import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

public class HpoOboParserTest2 {

	@Test
	public void parserTest() throws IOException {
		String input = ""; 
		PhenoMarksParser parser = new PhenoMarksParser();
		
		input = "mental retardation severe many x-linked several mild moderate developmental diffuse lesions hypotonia acute progressive dystrophy breathy failure unilateral hypoplastic underdeveloped skin aplasia (Syndactyly); (nystagmus) (cataract) (glaucoma) tralala";
		//input = "underdeveloped skin aplasia (Syndactyly); (nystagmus) (cataract) (glaucoma)";
		String result = parser.tagInput(input).getMarkedText();
		System.out.println(result);
		assertTrue(!result.contains("HP:0100754"));		
		assertTrue(result.contains("HP:0000518"));
		assertTrue(result.contains("HP:0000501"));
		
	}
	
}
